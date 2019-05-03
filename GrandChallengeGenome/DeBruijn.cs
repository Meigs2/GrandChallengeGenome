using GrandChallengeGenome.Models;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;

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
        public void AssembleGenome(int kMer)
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

            CleanGraph();

            ResolveContigs();

        }

        /// <summary>
        /// This function take the Graph we created and tries to remove paths that end early and paths that diverge.
        /// </summary>
        public void CleanGraph()
        {
            Console.WriteLine("Simplifying graph to improve memory and compute performance...");

            _contigGraph = _contigDictionaryGraph.Values.ToList();
            _contigDictionaryGraph = null;

            // Reduce pipes that start at a contig with one output ( a starting node )

            foreach (var contig in _contigGraph.Where(c => c.PreviousContigModels.Count == 0))
            {
                SimplifyForwardFromContig(contig);
            }
            DeleteMarked();

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

            // Reduce pipes that start after a merge.
            foreach (var afterConvergeContig in _contigGraph.Where(c => c.PreviousContigModels.Count > 1).ToList())
            {
                afterConvergeContig.Contig = afterConvergeContig.Contig[^1].ToString();
                SimplifyForwardFromContig(afterConvergeContig);
            }
            DeleteMarked();

            CountStartsMergesDiverges(_contigGraph);
        }

        private void ResolveContigs()
        {
            var localContigGraph = new List<ContigModel>(_contigGraph);
            var finalList = new List<ContigModel>();
            var startingContigs = localContigGraph
                .Where(c => c.PreviousContigModels.Count == 0).ToList();
            var routeList = new List<List<ContigModel>>();

            int currentStartingIndex = 0;
            // Locally Greedy DFS.
            foreach (var startingContig in startingContigs)
            {
                var currentRoute = new List<ContigModel>();
                var stack = new Stack<ContigModel>();
                stack.Push(startingContig);
                while (stack.Count != 0)
                {
                    var currentContig = stack.Pop();
                    currentRoute.Add(currentContig);

                    // If we've visited this contig already, keep searching.
                    if (currentContig.Visited)
                    {
                        var currentIndex = currentRoute.IndexOf(currentContig);
                        currentRoute.RemoveRange(currentIndex, currentRoute.Count - currentIndex);
                        continue;
                    }

                    // If we've reached an end contig
                    if (currentContig.NextContigModels.Count == 0)
                    {
                        // We've found a end contig! huzzah!

                        // If we've already visited this node, decide if we want to replace the route with our
                        // longer current route.
                        if (currentContig.VisitedEndNode)
                        {
                            var contestingRoute = routeList.FirstOrDefault(c => c.Contains(currentContig));

                            // get length of contesting route
                            var contigLength = 0;
                            foreach (var node in contestingRoute)
                            {
                                contigLength += node.Contig.Length;
                            }

                            // get length of our current route
                            var currentContigLength = 0;
                            foreach (var node in currentRoute)
                            {
                                currentContigLength += node.Contig.Length;
                            }

                            // if we're longer, replace
                            if (currentContigLength > contigLength)
                            {
                                routeList.Remove(contestingRoute);
                                routeList.Add(currentRoute);
                                break;
                            }

                            // else continue searching.
                            currentRoute.Remove(currentContig);
                            continue;
                        }

                        currentContig.VisitedEndNode = true;
                        foreach (var contigModel in currentRoute)
                        {
                            contigModel.IsPartOfRoute = true;
                        }
                        routeList.Add(currentRoute);
                        break;
                    }

                    // Resolve what to do if we've already been down this path.
                    if (currentContig.IsPartOfRoute)
                    {
                        var contestingRoute = routeList.FirstOrDefault(c => c.Contains(currentContig));

                        // Get length of contesting route up until this point.
                        var contigLength = 0;
                        foreach (var node in contestingRoute)
                        {
                            if (node == currentContig)
                            {
                                break;
                            }
                            contigLength += node.Contig.Length;
                        }

                        // Get length of our route up until this point.
                        var currentContigLength = 0;
                        foreach (var node in currentRoute)
                        {
                            if (node == currentContig)
                            {
                                break;
                            }
                            currentContigLength += node.Contig.Length;
                        }

                        // check which is longer, and replace if current is longer
                        if (currentContigLength > contigLength)
                        {
                            foreach (var contigModel in contestingRoute)
                            {
                                contigModel.IsPartOfRoute = false;
                            }

                            foreach (var contigModel in currentRoute)
                            {
                                contigModel.IsPartOfRoute = false;
                            }
                            var currentIndexInContesting = contestingRoute.IndexOf(currentContig);
                            contestingRoute.RemoveRange(0, currentIndexInContesting);
                            contestingRoute.InsertRange(0, currentRoute);
                            foreach (var contigModel in contestingRoute)
                            {
                                contigModel.IsPartOfRoute = true;
                            }
                            break;
                        }

                        // if not longer, continue searching elsewhere.
                        currentRoute.Remove(currentContig);
                        continue;
                    }

                    // Keep searching
                    currentContig.Visited = true;
                    foreach (var contigModel in currentContig.NextContigModels.OrderBy(x => x.Value))
                    {
                        stack.Push(contigModel.Key);
                    }
                }

                // When we decide to break, reset all the nodes that have been visited for the next search.
                // Clear visited nodes.
                var routeModels = localContigGraph.Where(c => c.IsPartOfRoute).ToList();
                foreach (var contigModel in routeModels)
                {
                    contigModel.Visited = false;
                }

                Console.WriteLine($"\rSearching starting point {currentStartingIndex++}/{startingContigs.Count}  ");
            }

            // Combine routes into single contigs.
            foreach (var route in routeList)
            {
                var combinedContig = new ContigModel();
                foreach (var contigModel in route)
                {
                    combinedContig.Contig += contigModel.Contig;
                }
                finalList.Add(combinedContig);
            }
            SaveData(finalList);
        }

        private void SaveData(List<ContigModel> contigsToExport)
        {
            List<string> exportStrings = new List<string>();
            foreach (var contigModel in contigsToExport)
            {
                exportStrings.Add($"> Built with k-mer of size {_kMerSize}.");
                exportStrings.Add(contigModel.Contig);
            }
            System.IO.File.WriteAllLines(Directory.GetCurrentDirectory() + $"\\rand.n.n50_{CalculateN50(contigsToExport)}.fa", exportStrings);
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
                bNextContigModel.PreviousContigModels.Add(a, b.NextContigModels[bNextContigModel]);
            }

            a.NextContigModels = b.NextContigModels;

            a.Contig = a.Contig + b.Contig[^1];
            b.PreviousContigModels = new Dictionary<ContigModel, int>();
            b.NextContigModels = new Dictionary<ContigModel, int>();
            b.Contig = string.Empty;
            b.MarkedForDeletion = true;
        }

        private void CountStartsMergesDiverges(List<ContigModel> contigs)
        {
            var baseContigs = new List<ContigModel>();
            var mergingContigs = new List<ContigModel>();
            var forkContigs = new List<ContigModel>();
            var endingContigs = new List<ContigModel>();
            foreach (var contig in contigs)
            {
                if (contig.PreviousContigModels.Count == 0) baseContigs.Add(contig);

                // diverging contigs
                if (contig.NextContigModels.Count > 1) forkContigs.Add(contig);

                // merging contigs
                if (contig.PreviousContigModels.Count > 1) mergingContigs.Add(contig);

                if (contig.NextContigModels.Count == 0)
                {
                    endingContigs.Add(contig);
                }
            }

            Console.WriteLine($"Found {baseContigs.Count} starting paths. Found {mergingContigs.Count} merging contigs. Found {forkContigs.Count} diverging contigs. Found {endingContigs.Count} end contigs.");
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
