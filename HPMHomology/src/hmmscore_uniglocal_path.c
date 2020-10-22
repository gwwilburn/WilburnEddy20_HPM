/* hpmscore_path.c */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_msa.h"
#include "esl_msafile.h"
#include "esl_sq.h"
#include "esl_vectorops.h"

#include "p7_config.h"
#include "hmmer.h"

static ESL_OPTIONS options[] = {
   /* name        type        default    env range togs  reqs  incomp            help                                            docgroup */
   { "-h",        eslARG_NONE,    FALSE, NULL, NULL, NULL, NULL, NULL,            "help; show brief info on version and usage",  1 },

   /* Options forcing which alphabet we're working in (normally autodetected) */
   { "--amino",   eslARG_NONE,    FALSE, NULL, NULL, NULL, NULL, "--dna,--rna",   "We are dealing with protein sequences",       2 },
   { "--rna",     eslARG_NONE,    FALSE, NULL, NULL, NULL, NULL, "--amino,--dna", "We are dealing with rna sequences",           2 },
   { "--dna",     eslARG_NONE,    FALSE, NULL, NULL, NULL, NULL, "--amino,--rna", "We are dealing with dna sequences",           2 },

   { 0,0,0,0,0,0,0,0,0,0 },
};

static char banner[] = "Score sequence(s) and path(s) uniglocally with an HMM";
static char usage[]  = "[-options] <hmmfile> <msafile> <score_outfile>";

static void
cmdline_failure(char *argv0, char *format, ...)
{
   va_list argp;
   printf("\nERROR: ");
   va_start(argp, format);
   vfprintf(stderr, format, argp);
   va_end(argp);
   esl_usage(stdout, argv0, usage);
   printf("\nTo see more help on available options, do %s -h\n\n", argv0);
   exit(1);
}

int main(int argc, char *argv[]){

   ESL_GETOPTS      *go        = NULL;
   ESL_ALPHABET     *abc       = NULL;
   char             *hmmfile   = NULL;                      /* input hmm filepath                        */
   char             *msafile   = NULL;                      /* input msa filepath                        */
   char             *scorefile = NULL;                      /* output score filepath                        */
   ESL_MSAFILE      *afp       = NULL;                      /* oopen input msa file stream               */
   P7_HMMFILE       *hfp       = NULL;                      /* open input hmm file stream     */
   ESL_MSA          *msa       = NULL;                      /* input msa                                 */
   P7_HMM           *hmm       = NULL;                      /* input hmm                      */
   ESL_SQ          **sq        = NULL;                      /* array of sequences                        */
   P7_TRACE        **tr        = NULL;
   P7_PROFILE       *gm        = NULL;                       /* H3 profile model                      */
   P7_BG            *bg        = NULL;                       /* H3 background model                   */
   FILE             *score_fp  = NULL;                      /* output score file parser object        */
   int              *matassign = NULL;                      /* array for ali match positions */
   float             hmmsc;                                 /* Joint hmm score of seq and path, in nats  */
   float             nullsc;                                /* Null transition score, in nats            */
   int               i;                                      /* msa position (column) index               */
   int               n;                                      /* msa sequence (row) index                  */
   int               fmt       = eslMSAFILE_STOCKHOLM;       /* input msa format #sorrynotsorry           */
   int               status;                                 /* easel return code                         */

   /* parse command line */
   go = esl_getopts_Create(options);
   if (esl_opt_ProcessCmdline(go, argc, argv) != eslOK) cmdline_failure(argv[0], "Failed to parse command line: %s\n", go->errbuf);
   if (esl_opt_VerifyConfig(go)               != eslOK) cmdline_failure(argv[0], "Error in app configuration:   %s\n", go->errbuf);

   /* output help if requested */
   if (esl_opt_GetBoolean(go, "-h") ) {
      esl_banner(stdout, argv[0], banner);
      esl_usage (stdout, argv[0], usage);
      puts("\n where options are:");
      esl_opt_DisplayHelp(stdout, go, 1, 2, 80);
      puts("\n Alphabet options:");
      esl_opt_DisplayHelp(stdout, go, 2, 2, 80);
      exit(0);
   }

   /* manually set alphabet if requested */
   if      (esl_opt_GetBoolean(go, "--amino")) abc = esl_alphabet_Create(eslAMINO);
   else if (esl_opt_GetBoolean(go, "--dna"))   abc = esl_alphabet_Create(eslDNA);
   else if (esl_opt_GetBoolean(go, "--rna"))   abc = esl_alphabet_Create(eslRNA);

   /* read arguments */
   if (esl_opt_ArgNumber(go) != 3)  cmdline_failure(argv[0], "Incorrect number of command line arguments.\n", go->errbuf);

   hmmfile   = esl_opt_GetArg(go, 1);
   msafile   = esl_opt_GetArg(go, 2);
   scorefile = esl_opt_GetArg(go, 3);

   /* Read in one HMM (and autodetect alphabet if not specified) */
   if (p7_hmmfile_OpenE(hmmfile, NULL, &hfp, NULL) != eslOK) esl_fatal("Failed to open HMM file %s", hmmfile);
   if (p7_hmmfile_Read(hfp, &abc, &hmm)            != eslOK) esl_fatal("Failed to read HMM");
   p7_hmmfile_Close(hfp);

   /* allow '-' and '.' as gap */
   esl_alphabet_SetEquiv(abc, '.', '-');

   /* read in MSA */
   if ((status = esl_msafile_Open(&abc, msafile, NULL, fmt, NULL, &afp)) != eslOK) esl_msafile_OpenFailure(afp, status);
   if ((status = esl_msafile_Read(afp, &msa))                            != eslOK) esl_fatal ("Failed to read MSA");
   esl_msafile_Close(afp);

   /* check to make sure MSA is digital and has RF line */
   if (! (msa->flags & eslMSA_DIGITAL)) esl_fatal("ERROR: need a digital msa");
   if (msa->rf == NULL)                 esl_fatal("ERROR: msa lacks an RF line");

   /* open score file */
   if ((score_fp = fopen(scorefile, "w")) == NULL) esl_fatal("Failed to open output hmm score file %s for writing", scorefile);

   /* set up profile and background models */
   bg = p7_bg_Create(abc);
   gm = p7_profile_Create(hmm->M, abc);

   /* configure profile in uniglocal mode */
   p7_profile_ConfigUniglocal(gm, hmm, bg, 400);

   /* extract sequences from MSA */
   ESL_ALLOC(sq, sizeof(ESL_SQ *) * msa->nseq);
   for (n = 0; n < msa->nseq; n++) {
      esl_sq_FetchFromMSA(msa, n, &sq[n]);
   }

   /* allocate memory for trace and match arrays */
   ESL_ALLOC(tr, sizeof(P7_TRACE *) * msa->nseq);
   ESL_ALLOC(matassign, sizeof(int) * (msa->alen + 1));
   esl_vec_ISet(matassign, msa->alen+1, 0);

   /* extract match states from alignment */
   for (i=0; i < msa->alen; i++){
      if (msa->rf[i] != '.') matassign[i+1] = 1;
      else                       matassign[i+1] = 0;
   }

   /* assign traces to msa seqs */
   p7_trace_FauxFromMSA(msa, matassign, p7_DEFAULT, tr);

   fprintf(score_fp, "id,logodds_path\n");

   /* loop over sequences and score */
   for (n = 0; n < msa->nseq; n++) {

      /* set profile and null models' target length */
      p7_bg_SetLength(bg, sq[n]->n);
      p7_profile_SetLength(gm, sq[n]->n);

      /* get partial hmm log-odds score, in nats */
      p7_trace_Score(tr[n], sq[n]->dsq, gm, &hmmsc);

      /* get null transition score, in nats */
      p7_bg_NullOne(bg, sq[n]->dsq, sq[n]->n, &nullsc);

      /* write full log-odds score, in bits, to output file */
      fprintf(score_fp, "%s,%.4f\n", sq[n]->name, (hmmsc - nullsc) / eslCONST_LOG2);
   }

   /* clean up */
   fclose(score_fp);
   for (n = 0; n < msa->nseq; n++) {
      p7_trace_Destroy(tr[n]);
      esl_sq_Destroy(sq[n]);
   }
   free(tr);
   free(sq);
   p7_hmm_Destroy(hmm);
   esl_msa_Destroy(msa);
   esl_alphabet_Destroy(abc);
   p7_profile_Destroy(gm);
   p7_bg_Destroy(bg);
   free(matassign);
   esl_getopts_Destroy(go);
   return 0;

   ERROR:
      return status;
}

