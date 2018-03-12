
//============================================================================//
//============================================================================//
//============================================================================//
//                                                                            //
//       CHANGE ONLY THE CODE INSIDE OF THE BODY OF "GenerateSubset".         //
//                                                                            //
//============================================================================//
//============================================================================//
//============================================================================//


#include <cstdio>
#include <cstdlib>
#include <vector>


// Holds the current minimum value found so far.
extern int global_min;

// Sum of all of the weights.
extern int total_sum;


//------------------------------------------------------------------------------
// GenerateSubset
//      Generates all possible subsets of weights recursively.
//
// PARAMETERS:
//      i       - Index of the element to be added or not added to the subset in
//                the current recursive step.
//      ans     - Set of weights in current subset.
//      S       - Reference to vector containing set of all weights.
//------------------------------------------------------------------------------
void GenerateSubset( int i, std::vector<int> ans, std::vector<int> &S ) {

//int sTotal = 0;
int ansTotal = 0;

    if ( i == S.size() ) {
//        for (int j =0;j<S.size();j++){sTotal+=S[j];}
        for (int j =0;j<ans.size();j++){ansTotal+=ans[j];} 
        int diff = abs(total_sum - 2*ansTotal);
        if (diff < global_min){
           #pragma omp critical
           if (diff < global_min){global_min=diff;}
        }

        // Execution reaches this point when a subset has been fully
        // constructed in "ans".
        // If this subset yields a better solution than had been previously
        // found, then update the result in "global_min".

        return;

    } else {

  #pragma omp task
        GenerateSubset(i+1, ans, S); // generate subsets without element i
  
  #pragma omp task
  {
        ans.push_back(S[i]); // add element i to the subset
  
        GenerateSubset(i+1, ans, S); // generate subsets with element i in them
  }



    }
}
