#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* easel includes */
#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_sq.h"
#include "esl_sqio.h"

/* h4 nwo incliudes */
#include "h4_hmmfile.h"
#include "h4_mode.h"
#include "h4_profile.h"
#include "h4_refmx.h"
#include "reference_dp.h"

/* declaration of internal functions */
int run_decoding(H4_PROFILE *hmm, ESL_SQ **sq, int nseq);

static ESL_OPTIONS options[] = {
   /* name          type            default env   range  togs  reqs            incomp           help                                                 docgroup */
   { "-h",          eslARG_NONE,    FALSE,  NULL, NULL,  NULL, NULL,           NULL,            "help; show brief info on version and usage",               1 },

   /* Options forcing which alphabet we're working in (normally autodetected) */
   { "--amino",     eslARG_NONE,    FALSE,  NULL, NULL,  NULL, NULL,           "--dna,--rna",   "<seqfile> contains protein sequences",                     2 },
   { "--rna",       eslARG_NONE,    FALSE,  NULL, NULL,  NULL, NULL,           "--dna,--amino", "<seqfile> contains RNA sequences",                         2 },
   { "--dna",       eslARG_NONE,    FALSE,  NULL, NULL,  NULL, NULL,           "--rna,--amino", "<seqfile> contains DNA sequences",                         2 },


   { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};

static char usage[]  = "[-options] <hmmfile> <seqfile>";
static char banner[] = "Examine hmm posterior decoding matrices";


static void
cmdline_failure(char *argv0, char *format, ...) {
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

   ESL_GETOPTS      *go            = NULL;                  /* application configuration                         */
   ESL_ALPHABET     *abc           = NULL;                  /* biological alphabet                               */
   char             *hmmfile       = NULL;                  /* input HMM filepath                                */
   char             *seqfile       = NULL;                  /* input seq filepath                                */
   H4_HMMFILE       *hmmfp         = NULL;                  /* open input hmm file stream                        */
   H4_PROFILE       *hmm           = NULL;                  /* hmm                                               */
   ESL_SQFILE       *sqfp          = NULL;                  /* open seq file stream                              */
   ESL_SQ          **sq            = NULL;                  /* array of sequences                                */
   int               format        = eslSQFILE_UNKNOWN;     /* input seq file format                             */
   int               n             = 0;                     /* sequence index                                    */
   int               nseq          = 0;                     /* number of seqs we deal with                       */
   int               status;                                /* Easel return code                                 */

   /* parse command line */
   go = esl_getopts_Create(options);
   if (esl_opt_ProcessCmdline(go, argc, argv) != eslOK) cmdline_failure(argv[0], "Failed to parse command line: %s\n", go->errbuf);
   if (esl_opt_VerifyConfig(go)               != eslOK) cmdline_failure(argv[0], "Error in app configuration:   %s\n", go->errbuf);

   if (esl_opt_GetBoolean(go, "-h") ) {
      esl_banner(stdout, argv[0], banner);
      esl_usage (stdout, argv[0], usage);
      puts("\n where options are:");
      esl_opt_DisplayHelp(stdout, go, 1, 2, 80);
      puts("\n Alphabet options:");
      esl_opt_DisplayHelp(stdout, go, 2, 2, 80);
      esl_getopts_Destroy(go);
      exit(0);
   }

   /* read arguments */
   if (esl_opt_ArgNumber(go) != 2)  cmdline_failure(argv[0], "Incorrect number of command line arguments.\n", go->errbuf);

   hmmfile   =  esl_opt_GetArg(go, 1);
   seqfile   =  esl_opt_GetArg(go, 2);

   /* if user has defined an alphabet we define it here */
   if        (esl_opt_GetBoolean(go, "--amino"))       abc = esl_alphabet_Create(eslAMINO);
   else if   (esl_opt_GetBoolean(go, "--rna"))         abc = esl_alphabet_Create(eslRNA);
   else if   (esl_opt_GetBoolean(go, "--dna"))         abc = esl_alphabet_Create(eslDNA);

   /* open the .hmm file */
   if ( h4_hmmfile_Open(hmmfile, NULL, &hmmfp) != eslOK) esl_fatal("couldn't open profile file %s", hmmfile);
   /* read first hmm  from hmm file */
   if ( h4_hmmfile_Read(hmmfp, &abc, &hmm)     != eslOK) esl_fatal("failed to read profile from file %s", hmmfile);
   h4_hmmfile_Close(hmmfp);

   /* open the sequence file */
   status = esl_sqfile_OpenDigital(abc, seqfile, format, NULL, &sqfp);
   if      (status == eslENOTFOUND) esl_fatal("No such file.");
   else if (status == eslEFORMAT)   esl_fatal("Format unrecognized.");
   else if (status == eslEINVAL)    esl_fatal("Can't autodetect stdin or .gz.");
   else if (status != eslOK)        esl_fatal("Open failed, code %d.", status);

   /* read sequences into array */
   ESL_REALLOC(sq, sizeof(ESL_SQ *));
   sq[n] = esl_sq_CreateDigital(abc);
   while ((status = esl_sqio_Read(sqfp, sq[n+nseq])) == eslOK) {
      nseq++;
      ESL_REALLOC(sq, sizeof(ESL_SQ *) * (nseq+1));
      sq[n+nseq] = esl_sq_CreateDigital(abc);
   }

   if (run_decoding(hmm, sq, nseq) != eslOK) esl_fatal("run_decoding() failed\n");


   fprintf(stdout, "hello world!\n");

   /* clean up and return */
   for (n = 0; n < nseq+1; n++) {
      esl_sq_Destroy(sq[n]);
   }
   free(sq);
   esl_sqfile_Close(sqfp);
   h4_profile_Destroy(hmm);
   esl_alphabet_Destroy(abc);
   esl_getopts_Destroy(go);
   return 0;

   ERROR:
      return status;

}

int run_decoding(H4_PROFILE *hmm, ESL_SQ **sq, int nseq) {

   H4_MODE        *mo      = h4_mode_Create();
   H4_REFMX       *fwd     = h4_refmx_Create(100, 100);
   H4_REFMX       *bck     = h4_refmx_Create(100, 100);
   H4_REFMX       *pp      = h4_refmx_Create(100, 100);
   float           fsc, bsc;
   int             n;
   int             i;
   int             k;
   float           mgk_sum, mgk_max, mgk_i;
   int             i_max;
   //float           i_sum;
   //int             status;
   char            errbuf[eslERRBUFSIZE];

   fprintf(stdout, "in run_decoding()\n");

   /* set HMM mode to uniglocal  */
   h4_mode_SetUniglocal(mo);

   for (n = 0; n < nseq; n++) {
      fprintf(stdout, "n: %d\n", n);

      h4_mode_SetLength(mo, sq[n]->n);
      h4_reference_Forward(sq[n]->dsq, sq[n]->n, hmm, mo, fwd, &fsc);
      //h4_refmx_Dump(stdout, fwd);
      printf("%s fwd %.6f\n", sq[n]->name, fsc);

      h4_reference_Backward(sq[n]->dsq, sq[n]->n, hmm, mo, bck, &bsc);
      printf("%s bck %.6f\n", sq[n]->name, bsc);

      h4_reference_Decoding(sq[n]->dsq, sq[n]->n, hmm, mo, fwd, bck, pp);

      if ( h4_refmx_Validate(pp, errbuf)          != eslOK) esl_fatal("%s\n", errbuf);

      // h4_refmx_Dump(stdout, pp);

      /* examine forward matrix */

      /* loop over global match states, find max and argmax */
      for (k = 1; k <= pp->M; k++) {

         mgk_sum = 0.;
         mgk_max = 0.;
         i_max = -1;
         for (i=1; i <= sq[n]->n; i++){
            mgk_i = pp->dp[i][k * h4R_NSCELLS + h4R_MG];

            if (mgk_i > mgk_max) {
               mgk_max = mgk_i;
               i_max = i;
            }

            mgk_sum += mgk_i;
         }
         fprintf(stdout, "k: %d, mgk_sum: %4f\n", k, mgk_sum);
      }

      h4_refmx_Reuse(fwd);
      h4_refmx_Reuse(bck);

      break;

   }



   /* clean up and return */
   h4_refmx_Destroy(pp);
   h4_refmx_Destroy(bck);
   h4_refmx_Destroy(fwd);
   h4_mode_Destroy(mo);

   return eslOK;
}
