using GrandChallengeGenome.Models;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.InteropServices;

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
        private List<ContigModel> _contigGraph = new List<ContigModel>();

        private int _kMerSize = 0;

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

            _baseContigModels = null;
        }

        /// <summary>
        /// Builds a DeBruijn Graph based on the contigs given in the constructor spaced by a specified k-mer.
        /// </summary>
        /// <param name="kMer"></param>
        public void BuildDeBruijnGraph(int kMer)
        {
            _kMerSize = kMer;
            Console.WriteLine($"Building Graph with k-mer of size {_kMerSize}...");

            // loop over contigs
            foreach (var contigModel in _contigModels)
            {
                ContigModel previousKMerContig = null;
                // split contig into k-mer size pieces
                for (var i = 0; i <= contigModel.Contig.Length - _kMerSize; i++)
                {
                    // Select range (C# 8.0 only, need to get the .net Core 3.0 Runtime and Development things to compile and use this)
                    var kMerContigSubstring = contigModel.Contig[i..i + _kMerSize];

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
            CountStartsMergesDiverges(_contigDictionaryGraph.Values.ToList());

            Console.WriteLine("Simplifying graph to improve memory and compute performance...");

            SimplifyGraph();

            CountStartsMergesDiverges(_contigGraph);

            //// combine contigs that are a "pipe" and have no other connections other than its neighbor and re-construct the graph.
            //// Doing this saves memory and compute performance down the road.
            //foreach (var startingContig in startingContigs)
            //{
            //    var currentContig = startingContig;
            //    while (currentContig.NextContigModels.Count != 0)
            //    {
            //        // If we have to investigate bubbles and splits.
            //        if (currentContig.NextContigModels.Count > 1)
            //        {

            //        }
            //        if (currentContig.PreviousContigModels.Count > 1)
            //        {

            //        }
            //        // If not, combine.
            //        else
            //        {

            //        }
            //    }
            //}
        }

        private void CountStartsMergesDiverges(List<ContigModel> contigs)
        {
            var baseContigs = new List<ContigModel>();
            var mergingContigs = new List<ContigModel>();
            var forkContigs = new List<ContigModel>();
            foreach (var contig in contigs)
            {
                if (contig.PreviousContigModels.Count == 0) baseContigs.Add(contig);

                // diverging contigs
                if (contig.NextContigModels.Count > 1) forkContigs.Add(contig);

                // merging contigs
                if (contig.PreviousContigModels.Count > 1) mergingContigs.Add(contig);
            }

            Console.WriteLine($"Found {baseContigs.Count} starting paths. Found {mergingContigs.Count} merging contigs. Found {forkContigs.Count} diverging contigs.");
        }

        /// <summary>
        /// This funtion simply finds all the kMer contigs with one in and out connection and makes them into one single contig.
        /// </summary>
        private void SimplifyGraph()
        {
            // Loop until we find a single in-out contig.
            _contigGraph = _contigDictionaryGraph.Values.ToList();
            _contigDictionaryGraph = null;

            // Begin process of merging pipe contigs

            // Reduce pipes that start at a contig with one output ( a starting node )
            var startingContigs = _contigGraph.Where(c => c.NextContigModels.Count == 1 && c.PreviousContigModels.Count == 0).ToList();
            foreach (var startingContig in startingContigs)
            {
                CondenceToNonPipe(startingContig);
            }
            startingContigs = null;
            _contigGraph.RemoveAll(c => c.MarkedForDeletion);

            // Reduce pipes that start at a convergence
            var convergingContigs = _contigGraph.Where(c => c.PreviousContigModels.Count > 1);
            foreach (var convergingContig in convergingContigs)
            {
                CondenceToNonPipe(convergingContig);
            }
            convergingContigs = null;
            _contigGraph.RemoveAll(c => c.MarkedForDeletion);

            // Reduce paths that start after a divergence (bubbles/branches)

            var divergingContigs = _contigGraph.Where(c => c.NextContigModels.Count > 1);
            foreach (var divergingContig in divergingContigs)
            {
                foreach (var forkedContig in divergingContig.NextContigModels)
                {
                    CondenceToNonPipe(forkedContig);
                }
            }
            divergingContigs = null;
            _contigGraph.RemoveAll(c => c.MarkedForDeletion);

            foreach (var contigModel in _contigGraph)
            {
                Console.WriteLine($"Found a contig with length {contigModel.Contig.Length}.");
            }
                
        }

        /// <summary>
        /// Search Given Pipe
        /// </summary>
        /// <param name="currentModel"></param>
        /// <returns></returns>
        private void CondenceToNonPipe(ContigModel currentModel)
        {
            while (true)
            {
                // if we have one connection, merge with next node, continue.
                if (currentModel.NextContigModels.Count == 1)
                {
                    // check if the next node has more than one previous contigs coming into it.
                    var nextModel = currentModel.NextContigModels.First();
                    if (nextModel.PreviousContigModels.Count > 1)
                    {
                        break;
                    }
                    MergeNodes(currentModel, nextModel);
                    continue;
                }
                break;
            }
        }

        // Merges A with b.
        private void MergeNodes(ContigModel a, ContigModel b)
        {
            a.Contig = a.Contig + b.Contig[^1];
            a.NextContigModels = b.NextContigModels;
            b.MarkedForDeletion = true;
        }

        public void DisposeMemory()
        {
            _baseContigModels = null;
            _contigDictionaryGraph = null;
            _contigModels = null;
            _contigGraph = null;
        }
    }
}
