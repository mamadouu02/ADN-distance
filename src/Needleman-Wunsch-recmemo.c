/**
 * \file Needleman-Wunsch-recmemo.c
 * \brief recursive implementation with memoization of Needleman-Wunsch global alignment algorithm that computes the distance between two genetic sequences 
 * \version 0.1
 * \date 03/10/2022 
 * \author Jean-Louis Roch (Ensimag, Grenoble-INP - University Grenoble-Alpes) jean-louis.roch@grenoble-inp.fr
 *
 * Documentation: see Needleman-Wunsch-recmemo.h
 * Costs of basic base opertaions (SUBSTITUTION_COST, SUBSTITUTION_UNKNOWN_COST, INSERTION_COST) are
 * defined in Needleman-Wunsch-recmemo.h
 */


#include "Needleman-Wunsch-recmemo.h"
#include <stdio.h>  
#include <stdlib.h> 
#include <string.h> /* for strchr */
// #include <ctype.h> /* for toupper */

#include "characters_to_base.h" /* mapping from char to base */

/*****************************************************************************/
   
/* Context of the memoization : passed to all recursive calls */
/** \def NOT_YET_COMPUTED
 * \brief default value for memoization of minimal distance (defined as an impossible value for a distance, -1).
 */
#define NOT_YET_COMPUTED -1L 

/** \def S
 * \brief threshold for cache-oblivious Needleman-Wunsch algorithm.
 */
#define S 16

/** \struct NW_MemoContext
 * \brief data for memoization of recursive Needleman-Wunsch algorithm 
*/
struct NW_MemoContext
{
    char *X ; /*!< the longest genetic sequences */
    char *Y ; /*!< the shortest genetic sequences */
    size_t M; /*!< length of X */
    size_t N; /*!< length of Y,  N <= M */
    long **memo; /*!< memoization table to store memo[0..M][0..N] (including stopping conditions phi(M,j) and phi(i,N) */
};

/** \struct NW_MemoContextIter
 * \brief data for memoization of iterative Needleman-Wunsch algorithm 
*/
struct NW_MemoContextIter
{
    char *X ; /*!< the longest genetic sequences */
    char *Y ; /*!< the shortest genetic sequences */
    size_t M; /*!< length of X */
    size_t N; /*!< length of Y,  N <= M */
    long *memo; /*!< memoization table to store memo[0..N] (including stopping conditions phi(i,N) */
};

/** \struct NW_MemoContextCO
 * \brief data for memoization of cache-oblivious Needleman-Wunsch algorithm 
*/
struct NW_MemoContextCO
{
    char *X ; /*!< the longest genetic sequences */
    char *Y ; /*!< the shortest genetic sequences */
    size_t M; /*!< length of X */
    size_t N; /*!< length of Y,  N <= M */
    long *rmemo; /*!< memoization table to store memo[0..N] (including stopping conditions phi(i,N) */
    long *cmemo; /*!< memoization table to store memo[0..S] */
};

/**
 *  static long EditDistance_NW_RecMemo(struct NW_MemoContext *c, size_t i, size_t j) 
 * \brief  EditDistance_NW_RecMemo :  Private (static)  recursive function with memoization \
 * direct implementation of Needleman-Wursch extended to manage FASTA sequences (cf TP description)
 * \param c : data passed for recursive calls that includes the memoization array 
 * \param i : starting position of the left sequence :  c->X[ i .. c->M ] 
 * \param j : starting position of the right sequence :  c->Y[ j .. c->N ] 
 */ 
static long EditDistance_NW_RecMemo(struct NW_MemoContext *c, size_t i, size_t j) 
/* compute and returns phi(i,j) using data in c -allocated and initialized by EditDistance_NW_Rec */
{
   if (c->memo[i][j] == NOT_YET_COMPUTED)
   {  
      long res ;
      char Xi = c->X[i] ;
      char Yj = c->Y[j] ;
      if (i == c->M) /* Reach end of X */
      {  if (j == c->N) res = 0;  /* Reach end of Y too */
         else res = (isBase(Yj) ? INSERTION_COST : 0) + EditDistance_NW_RecMemo(c, i, j+1) ;
      }
      else if (j == c->N) /* Reach end of Y but not end of X */
      {  res = (isBase(Xi) ? INSERTION_COST : 0) + EditDistance_NW_RecMemo(c, i+1, j) ;
      }
      else if (! isBase(Xi))  /* skip ccharacter in Xi that is not a base */
      {  ManageBaseError( Xi ) ;
         res = EditDistance_NW_RecMemo(c, i+1, j) ;
      } 
      else if (! isBase(Yj))  /* skip ccharacter in Yj that is not a base */
      {  ManageBaseError( Yj ) ;
         res = EditDistance_NW_RecMemo(c, i, j+1) ;
      }
      else  
      {  /* Note that stopping conditions (i==M) and (j==N) are already stored in c->memo (cf EditDistance_NW_Rec) */ 
         long min = /* initialization  with cas 1*/
                   ( isUnknownBase(Xi) ?  SUBSTITUTION_UNKNOWN_COST 
                          : ( isSameBase(Xi, Yj) ? 0 : SUBSTITUTION_COST ) 
                   )
                   + EditDistance_NW_RecMemo(c, i+1, j+1) ; 
         { long cas2 = INSERTION_COST + EditDistance_NW_RecMemo(c, i+1, j) ;      
           if (cas2 < min) min = cas2 ;
         }
         { long cas3 = INSERTION_COST + EditDistance_NW_RecMemo(c, i, j+1) ;      
           if (cas3 < min) min = cas3 ; 
         }
         res = min ;
      }
       c->memo[i][j] = res ;
   }
   return c->memo[i][j] ;
}

/* EditDistance_NW_Rec :  is the main function to call, cf .h for specification 
 * It allocates and initailizes data (NW_MemoContext) for memoization and call the 
 * recursivefunction EditDistance_NW_RecMemo 
 * See .h file for documentation
 */
long EditDistance_NW_Rec(char* A, size_t lengthA, char* B, size_t lengthB)
{
   _init_base_match() ;
   struct NW_MemoContext ctx;
   if (lengthA >= lengthB) /* X is the longest sequence, Y the shortest */
   {  ctx.X = A ;
      ctx.M = lengthA ;
      ctx.Y = B ;
      ctx.N = lengthB ;
   }
   else
   {  ctx.X = B ;
      ctx.M = lengthB ;
      ctx.Y = A ;
      ctx.N = lengthA ;
   }
   size_t M = ctx.M ;
   size_t N = ctx.N ;
   {  /* Allocation and initialization of ctx.memo to NOT_YET_COMPUTED*/
      /* Note: memo is of size (N+1)*(M+1) but is stored as (M+1) distinct arrays each with (N+1) continuous elements 
       * It would have been possible to allocate only one big array memezone of (M+1)*(N+1) elements 
       * and then memo as an array of (M+1) pointers, the memo[i] being the address of memzone[i*(N+1)].
       */ 
      ctx.memo = (long **) malloc ( (M+1) * sizeof(long *)) ;
      if (ctx.memo == NULL) { perror("EditDistance_NW_Rec: malloc of ctx_memo" ); exit(EXIT_FAILURE); }
      for (int i=0; i <= M; ++i) 
      {  ctx.memo[i] = (long*) malloc( (N+1) * sizeof(long));
         if (ctx.memo[i] == NULL) { perror("EditDistance_NW_Rec: malloc of ctx_memo[i]" ); exit(EXIT_FAILURE); }
         for (int j=0; j<=N; ++j) ctx.memo[i][j] = NOT_YET_COMPUTED ;
      }
   }    
   
   /* Compute phi(0,0) = ctx.memo[0][0] by calling the recursive function EditDistance_NW_RecMemo */
   long res = EditDistance_NW_RecMemo( &ctx, 0, 0 ) ;
    
   { /* Deallocation of ctx.memo */
      for (int i=0; i <= M; ++i) free( ctx.memo[i] ) ;
      free( ctx.memo ) ;
   }
   return res ;
}

long EditDistance_NW_Iter(char *A, size_t lengthA, char *B, size_t lengthB)
{
   _init_base_match();
   struct NW_MemoContextIter ctx;

   if (lengthA >= lengthB) {
      ctx.X = A;
      ctx.M = lengthA;
      ctx.Y = B;
      ctx.N = lengthB;
   } else {
      ctx.X = B;
      ctx.M = lengthB;
      ctx.Y = A;
      ctx.N = lengthA;
   }

   size_t M = ctx.M;
   size_t N = ctx.N;

   ctx.memo = malloc(sizeof(long[N+1]));

   if (ctx.memo == NULL) {
      perror("EditDistance_NW_Iter: malloc of ctx.memo");
      exit(EXIT_FAILURE);
   }

   ctx.memo[N] = 0;

   for (int j = N-1; j >= 0; --j) {
      ctx.memo[j] = (isBase(ctx.Y[j]) ? INSERTION_COST : 0) + ctx.memo[j+1];
   }

   for (int i = M-1; i >= 0; --i) {
      long tmp = ctx.memo[N];
      ctx.memo[N] = (isBase(ctx.X[i]) ? INSERTION_COST : 0) + tmp;

      for (int j = N-1; j >= 0; --j) {
         if (!isBase(ctx.X[i])) {
            ManageBaseError(ctx.X[i]);
            tmp = ctx.memo[j];
         } else if (!isBase(ctx.Y[j])) {
            ManageBaseError(ctx.Y[j]);
            tmp = ctx.memo[j];
            ctx.memo[j] = ctx.memo[j+1];
         } else {
            long cost1 = sigma(ctx.X[i], ctx.Y[j]) + tmp;
            long cost2 = INSERTION_COST + ctx.memo[j];
            long cost3 = INSERTION_COST + ctx.memo[j+1];
            tmp = ctx.memo[j];
            ctx.memo[j] = MIN(cost1, cost2, cost3);
         }
      }
   }

   long res = ctx.memo[0];
   free(ctx.memo);

   return res;
}

 long EditDistance_NW_CA(char *A, size_t lengthA, char *B, size_t lengthB, int Z)
{
   _init_base_match();
   struct NW_MemoContextIter ctx;

   if (lengthA >= lengthB) {
      ctx.X = A;
      ctx.M = lengthA;
      ctx.Y = B;
      ctx.N = lengthB;
   } else {
      ctx.X = B;
      ctx.M = lengthB;
      ctx.Y = A;
      ctx.N = lengthA;
   }

   size_t M = ctx.M;
   size_t N = ctx.N;

   ctx.memo = malloc(sizeof(long[N+1]));

   if (ctx.memo == NULL) {
      perror("EditDistance_NW_CA: malloc of ctx.memo");
      exit(EXIT_FAILURE);
   }

   ctx.memo[N] = 0;
   
   for (int j = N-1; j >= 0; --j) {
      ctx.memo[j] = (isBase(ctx.Y[j]) ? INSERTION_COST : 0) + ctx.memo[j+1];
   }

   int K = Z / 64;

   for (int I = M-1; I >= 0; I -= K) {
      int i_end = (0 >= I - K) ? 0 : I - K;

      for (int J = N-1; J >= 0; J -= K) {
         int j_end = (0 >= J - K) ? 0: J - K;

         for (int i = I; i >= i_end; --i) {
            long tmp = ctx.memo[N];
            ctx.memo[N] = (isBase(ctx.X[i]) ? INSERTION_COST : 0) + tmp;

            for (int j = J; j >= j_end; --j) {
               if (!isBase(ctx.X[i])) {
                  ManageBaseError(ctx.X[i]);
                  tmp = ctx.memo[j];
               } else if (!isBase(ctx.Y[j])) {
                  ManageBaseError(ctx.Y[j]);
                  tmp = ctx.memo[j];
                  ctx.memo[j] = ctx.memo[j+1];
               } else {
                  long cost1 = sigma(ctx.X[i], ctx.Y[j]) + tmp;
                  long cost2 = INSERTION_COST + ctx.memo[j];
                  long cost3 = INSERTION_COST + ctx.memo[j+1];
                  tmp = ctx.memo[j];
                  ctx.memo[j] = MIN(cost1, cost2, cost3);
               }
            }
         }
      }
   }

   long res = ctx.memo[0];
   free(ctx.memo);

   return res;
}

void blocking(struct NW_MemoContextCO *c, int i_start, int i_stop)
{
   if (i_start - i_stop < S) {
      for (int i = i_start; i >= i_stop; --i) {
         int k = i - i_stop;

         if (i == c->M) {
            c->cmemo[k] = 0;
         } else {
            c->cmemo[k] = (isBase(c->X[i]) ? INSERTION_COST : 0) + ((i == i_start) ? c->rmemo[c->N] : c->cmemo[k + 1]);
         }
      }

      for (int j = c->N-1; j >= 0; --j) {
         for (int i = i_start; i >= i_stop; --i) {
            int k = i - i_stop;

            if (i == c->M) {
               c->rmemo[j+1] = c->cmemo[k];
               c->cmemo[k] = (isBase(c->Y[j]) ? INSERTION_COST : 0) + c->cmemo[k + 1];
            } else {
               if (!isBase(c->X[i])) {
                  ManageBaseError(c->X[i]);
                  c->rmemo[j+1] = c->cmemo[k];
                  c->cmemo[k] = ((i == i_start) ? c->rmemo[j] : c->cmemo[k + 1]);
               } else if (!isBase(c->Y[j])) {
                  ManageBaseError(c->Y[j]);
                  c->rmemo[j+1] = c->cmemo[k];
               } else {
                  long cost1 = sigma(c->X[i], c->Y[j]) + c->rmemo[j+1];
                  long cost2 = INSERTION_COST + ((i == i_start) ? c->rmemo[j] : c->cmemo[k + 1]);
                  long cost3 = INSERTION_COST + c->cmemo[k];
                  c->rmemo[j+1] = c->cmemo[k];
                  c->cmemo[k] = MIN(cost1, cost2, cost3);
               }
            }
         }
      }

      c->rmemo[0] = c->cmemo[0];
   } else {
      int i_mid = (i_start + i_stop) / 2;
      blocking(c, i_start, i_mid + 1);
      blocking(c, i_mid, i_stop);
   }
}

long EditDistance_NW_CO(char *A, size_t lengthA, char *B, size_t lengthB)
{
   _init_base_match();
   struct NW_MemoContextCO ctx;

   if (lengthA >= lengthB) {
      ctx.X = A;
      ctx.M = lengthA;
      ctx.Y = B;
      ctx.N = lengthB;
   } else {
      ctx.X = B;
      ctx.M = lengthB;
      ctx.Y = A;
      ctx.N = lengthA;
   }

   size_t M = ctx.M;
   size_t N = ctx.N;

   ctx.rmemo = malloc(sizeof(long[N+1]));

   if (ctx.rmemo == NULL) {
      perror("EditDistance_NW_Iter: malloc of ctx.rmemo");
      exit(EXIT_FAILURE);
   }

   ctx.cmemo = malloc(sizeof(long[S]));

   if (ctx.cmemo == NULL) {
      perror("EditDistance_NW_Iter: malloc of ctx.cmemo");
      exit(EXIT_FAILURE);
   }

   blocking(&ctx, M, 0);

   long res = ctx.rmemo[0];
   free(ctx.rmemo);
   free(ctx.cmemo);

   return res;
}
