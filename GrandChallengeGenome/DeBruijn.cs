using GrandChallengeGenome.Models;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.InteropServices;
using Microsoft.VisualBasic.FileIO;

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

        private List<ContigModel> _finalContigs = new List<ContigModel>();

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
                    var kMerContigSubstring = contigModel.Contig[i..i + _kMerSize - 1];

                    // Check if our kMerContig is in the dictionary already
                    _contigDictionaryGraph.TryGetValue(kMerContigSubstring, out var currentKMerContig);

                    // add contig to dict if not already included, and re-set current contig to make sure current contig is the one in the dictionary
                    if (currentKMerContig == null)
                    {
                        _contigDictionaryGraph.Add(kMerContigSubstring, new ContigModel() { Contig = kMerContigSubstring });
                        currentKMerContig = _contigDictionaryGraph[kMerContigSubstring];
                    }

                    // if we're looking at the first kmer substring of larger contig, set previous to current and go to the next substring.
                    if (previousKMerContig == null)
                    {
                        previousKMerContig = currentKMerContig;
                        continue;
                    }

                    if (previousKMerContig.NextContigModels.TryGetValue(currentKMerContig, out var nextValue))
                    {
                        previousKMerContig.NextContigModels[currentKMerContig] += 1;
                        currentKMerContig.PreviousContigModels[previousKMerContig] += 1;
                    }
                    else
                    {
                        previousKMerContig.NextContigModels.Add(currentKMerContig, 1);
                        currentKMerContig.PreviousContigModels.Add(previousKMerContig, 1);
                    }

                    previousKMerContig = currentKMerContig;
                }
            }

            Console.WriteLine("Graph Built Successfully!");

            CleanupGraph();
        }

        /// <summary>
        /// This function take the Graph we created and tries to remove paths that end early and paths that diverge.
        /// </summary>
        public void CleanupGraph()
        {
            Console.WriteLine("Simplifying graph to improve memory and compute performance...");

            SimplifyGraph();
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
            _contigGraph = _contigDictionaryGraph.Values.ToList();
            _contigDictionaryGraph = null;

            ResolveBubblesAndBranches();

        }

        private void ResolveBubblesAndBranches()
        {
            CountStartsMergesDiverges(_contigGraph);

            // Reduce pipes that start at a contig with one output ( a starting node )

            foreach (var contig in _contigGraph.Where(c => c.PreviousContigModels.Count == 0))
            {
                SimplifyForwardFromContig(contig);
            }
            DeleteMarked();
            //startingContigs = null;


            var startingContigs = _contigGraph
            .Where(c => c.NextContigModels.Count == 1 && c.PreviousContigModels.Count == 0).ToList();

            CountStartsMergesDiverges(_contigGraph);

            // Reduce pipes that start after a diverge
            foreach (var afterDivergeContig in _contigGraph.Where(c => c.NextContigModels.Count > 1))
            {
                foreach (var contig in afterDivergeContig.NextContigModels.Keys)
                {
                    contig.Contig = contig.Contig[^1].ToString();
                    SimplifyForwardFromContig(contig);
                }
            }
            DeleteMarked();
            CountStartsMergesDiverges(_contigGraph);



            // Reduce pipes that start after a merge.

            foreach (var afterConvergeContig in _contigGraph.Where(c => c.PreviousContigModels.Count > 1).ToList())
            {
                afterConvergeContig.Contig = afterConvergeContig.Contig[^1].ToString();
                SimplifyForwardFromContig(afterConvergeContig);
            }
            DeleteMarked();

            //var a = SearchPathsForwardGreedy();

            // Reverse search from ending nodes to beginning ones, until we reach a starting node.
            //afterConvergeContigs = null;

            //afterConvergeContigs = null;
        }

        private void SearchPathsForwardGreedy()
        {

        }

        private void DeleteMarked()
        {
            _contigGraph.RemoveAll(c => c.MarkedForDeletion);
        }

        /// <summary>
        /// Recursively simplify connections.
        /// </summary>
        /// <param name="currentModel"></param>
        /// <returns></returns>
        private void SimplifyForwardFromContig(ContigModel currentModel)
        {
            while (true)
            {
                // if we currently have at least one connection
                if (currentModel.NextContigModels.Count >= 1)
                {
                    var nextModel = currentModel.NextContigModels.Keys.First();

                    // Stop if the next contig has more than once connection to it.
                    if (nextModel.PreviousContigModels.Count > 1)
                    {
                        break;
                    }

                    // Stop if we reach a split in the graph.
                    if (currentModel.NextContigModels.Count > 1)
                    {
                        break;
                    }

                    // We merge if next contigs is part of the pipe.
                    MergeNodes(currentModel, nextModel);

                    continue;
                }
                // if we reach the end of a contig, stop
                break;
            }
        }

        // Merges A with b.
        private void MergeNodes(ContigModel a, ContigModel b)
        {
            // carry over visit count
            foreach (var bNextContigModel in b.NextContigModels.Keys)
            {
                var success = bNextContigModel.PreviousContigModels.Remove(b);
                bNextContigModel.PreviousContigModels.Add(a,b.NextContigModels[bNextContigModel]);
            }

            a.NextContigModels = b.NextContigModels;

            a.Contig = a.Contig + b.Contig[^1];
            b.PreviousContigModels = new Dictionary<ContigModel, int>();
            b.NextContigModels = new Dictionary<ContigModel, int>();
            b.Contig = string.Empty;
            b.MarkedForDeletion = true;
        }

        public void DisposeMemory()
        {
            _baseContigModels = null;
            _contigDictionaryGraph = null;
            _contigModels = null;
            _contigGraph = null;
        }
        public int CalculateN50(List<ContigModel> contigsList)
        {
            contigsList = contigsList.OrderByDescending(c => c.Contig.Length).ToList();
            var totalLength = 0;
            foreach (var contigModel in contigsList)
            {
                totalLength += contigModel.Contig.Length;
            }
            var runningSum = 0;
            foreach (var contigModel in contigsList)
            {
                runningSum = runningSum += contigModel.Contig.Length;
                if ((double)runningSum >= (double)totalLength / 2)
                {
                    return contigModel.Contig.Length;
                }
            }

            return -1;
        }
    }
}
