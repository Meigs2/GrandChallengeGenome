using GrandChallengeGenome.Models;
using System;
using System.Collections.Generic;

namespace GrandChallengeGenome
{
    /// <summary>
    /// This class does everything related to De Bruijin graphs.
    /// </summary>
    public class DeBruijn
    {
        private List<BaseContigModel> _baseContigModels = new List<BaseContigModel>();
        private List<ContigModel> _contigModels = new List<ContigModel>();

        /// <summary>
        /// This dictionary stores the kMer sized substring and returns the related contig
        /// </summary>
        private Dictionary<string, ContigModel> _contigDictionaryGraph = new Dictionary<string, ContigModel>();

        public DeBruijn(List<BaseContigModel> baseContigModels)
        {
            _baseContigModels = baseContigModels;

            // Convert our base models to normal ContigModels
            ConvertBaseModelsToNormalModels();
        }

        private void ConvertBaseModelsToNormalModels()
        {
            // Convert our base models
            foreach (var baseContigModel in _baseContigModels)
            {
                var model = new ContigModel();
                model.Contig = baseContigModel.Contig;
                _contigModels.Add(model);
            }
        }

        /// <summary>
        /// Builds a DeBruijn Graph based on the contigs given in the constructor spaced by a specified k-mer.
        /// </summary>
        /// <param name="kMer"></param>
        public void BuildDeBruijnGraph(int kMer)
        {
            Console.WriteLine($"Building Graph with k-mer of size {kMer}...");

            // loop over contigs
            foreach (var contigModel in _contigModels)
            {
                ContigModel previousKMerContig = null;
                // split contig into k-mer size pieces
                for (var i = 0; i < contigModel.Contig.Length - kMer; i++)
                {
                    // Select range (C# 8.0 only, need to get the .net Core 3.0 Runtime and Development things to compile and use this)
                    var kMerContigSubstring = contigModel.Contig[i..i + kMer];

                    // Check if our kMerContig is in the dictionary already
                    _contigDictionaryGraph.TryGetValue(kMerContigSubstring, out var currentKMerContig);

                    // add contig to dict if not already included, and re-set current contig to make sure current contig is the one in the dictionary
                    if (currentKMerContig == null)
                    {
                        _contigDictionaryGraph.Add(kMerContigSubstring, new ContigModel(){Contig = kMerContigSubstring});
                        currentKMerContig = _contigDictionaryGraph[kMerContigSubstring];
                    }

                    // if we're looking at the first kmer substring of larger contig, set previous to current and go to the next substring.
                    if (previousKMerContig == null)
                    {
                        previousKMerContig = currentKMerContig;
                        continue;
                    }

                    // link current contig to previous contig, and vice-versa
                    currentKMerContig.PreviousContigModels.Add(previousKMerContig);
                    previousKMerContig.NextContigModels.Add(currentKMerContig);

                    previousKMerContig = currentKMerContig;
                }
            }

            Console.WriteLine("Graph Build Successfully!");
            Console.WriteLine("Attempting to clean up the newly created graph...");
            CleanupGraph();
        }

        /// <summary>
        /// This function take the Graph we created and tries to remove paths that end early and paths that diverge.
        /// </summary>
        public void CleanupGraph()
        {
            var startingContigs = new List<ContigModel>();
            foreach (var contig in _contigDictionaryGraph)
            {
                // traverse graph starting where there are NO 
                if (contig.Value.PreviousContigModels.Count == 0)
                {
                    startingContigs.Add(contig.Value);
                }
            }

            Console.WriteLine($"Found {startingContigs.Count} starting paths. Combining paths to improve memory and compute performance...");

            // combine contigs that are a "pipe" and have no other connections other than its neighbor and re-construct the graph.
            // Doing this saves memory and compute performance down the road.
            foreach (var startingContig in startingContigs)
            {
                var currentContig = startingContig;
                while (currentContig.NextContigModels.Count != 0)
                {
                    // If we have to investigate bubbles and splits.
                    if (currentContig.NextContigModels.Count > 1)
                    {
                        
                    }
                    if (currentContig.PreviousContigModels.Count > 1)
                    {
                        
                    }
                    // If not, combine.
                    else
                    {
                        
                    }
                }
            }
        }

        private void SearchSplit(ContigModel currentModel)
        {

        }
    }
}
