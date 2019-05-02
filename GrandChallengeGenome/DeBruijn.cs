using GrandChallengeGenome.Models;
using System;
using System.Collections.Generic;
using System.ComponentModel.Design;
using System.IO;
using System.Linq;
using System.Net.Mime;
using System.Reflection;
using System.Runtime.InteropServices;
using Microsoft.VisualBasic.FileIO;
using SearchOption = Microsoft.VisualBasic.FileIO.SearchOption;

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

            GC.Collect();
            GC.WaitForPendingFinalizers();

            ResolveBubblesAndBranches();

            ProduceFinalContigs();

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

            // Reduce num starting nodes

            // Reduce num ending nodes
        }

        private void ProduceFinalContigs()
        {
            var localContigGraph = new List<ContigModel>(_contigGraph);

            var finalList = new List<ContigModel>();

            var startingContigs = localContigGraph
                .Where(c => c.PreviousContigModels.Count == 0).ToList();

            //ContigModel bestContig;

            var routeList = new List<List<ContigModel>>();

            // DFS for 
            foreach (var startingContig in startingContigs)
            {
                // Recurtsive DFS on node
                var route = new List<ContigModel>();
                var stack = new Stack<ContigModel>();
                stack.Push(startingContig);
                while (stack.Count != 0)
                {
                    var currentContig = stack.Pop();
                    route.Add(currentContig);
                    if (currentContig.Visited)
                    {
                        var currentIndex = route.IndexOf(currentContig);
                        route.RemoveRange(currentIndex, route.Count - currentIndex);
                        continue;
                    }

                    if (currentContig.NextContigModels.Count == 0)
                    {
                        // We've found a end node! huzzah!

                        if (currentContig.VisitedEndNode)
                        {
                            var contestingRoute = routeList.FirstOrDefault(c => c.Contains(currentContig));
                            // get length of consesting route
                            var contestingString = "";
                            foreach (var node in contestingRoute)
                            {
                                contestingString += node.Contig;
                            }

                            // get length of our current route
                            var currentString = "";
                            foreach (var node in route)
                            {
                                currentString += node.Contig;
                            }

                            if (currentString.Length > contestingString.Length)
                            {
                                routeList.Remove(contestingRoute);
                                routeList.Add(route);
                                break;
                            }

                            // check which one is longer, and replace, etc.
                            route.Remove(currentContig);
                            continue;
                        }

                        currentContig.VisitedEndNode = true;
                        routeList.Add(route);
                        break;
                    }

                    currentContig.Visited = true;
                    var test = currentContig.NextContigModels.OrderBy(x => x.Value).ToList();
                    foreach (var contigModel in currentContig.NextContigModels.OrderBy(x => x.Value))
                    {
                        stack.Push(contigModel.Key);
                    }
                }
                // Clear visited nodes.
                var a = localContigGraph.Where(c => c.Visited);
                foreach (var contigModel in a)
                {
                    contigModel.Visited = false;
                }
            }

            foreach (var route in routeList)
            {
                var combinedContig = new ContigModel(); 
                foreach (var contigModel in route)
                {
                    combinedContig.Contig += contigModel.Contig;
                }
                finalList.Add(combinedContig);
            }

            //// Lets try a greedy approach to finding the shortest path, then checking if there are multiple instances of the
            //// Same graph in the new list.
            //foreach (var startingContig in startingContigs)
            //{
            //    bestContig = startingContig;
            //    var currentContig = startingContig;
            //    while (currentContig.NextContigModels.Count != 0)
            //    {
            //        // if split, pick path with most visits
            //        if (currentContig.NextContigModels.Count > 1)
            //        {
            //            KeyValuePair<ContigModel, int> mostVisitedContig = new KeyValuePair<ContigModel, int>();
            //            foreach (var nextContig in currentContig.NextContigModels)
            //            {
            //                if (mostVisitedContig.Value == 0 || nextContig.Value > mostVisitedContig.Value)
            //                {
            //                    if (nextContig.Key.Visited)
            //                    {
            //                        goto exitIfVisited;
            //                    }

            //                    mostVisitedContig = nextContig;
            //                }
            //            }

            //            bestContig.Contig += mostVisitedContig.Key.Contig;
            //            mostVisitedContig.Key.Visited = true;
            //            currentContig = mostVisitedContig.Key;
            //            continue;

            //            exitIfVisited:
            //            bestContig.Contig += mostVisitedContig.Key.Contig;
            //            break;

            //        }
            //        if (currentContig.NextContigModels.Keys.First().Visited == false)
            //        {
            //            bestContig.Contig += currentContig.NextContigModels.Keys.First().Contig;
            //            currentContig.NextContigModels.Keys.First().Visited = true;
            //            currentContig = currentContig.NextContigModels.Keys.First();
            //        }
            //        else
            //        {
            //            break;
            //        }
            //    }


            //    finalList.Add(bestContig);
            //    // Reset visited Graph.
            //    var a = localContigGraph.Where(c => c.Visited);
            //    foreach (var contigModel in a)
            //    {
            //        contigModel.Visited = false;
            //    }
            //}

            //var copy = finalList.ToList();
            //finalList.RemoveAll(x => copy.Any(y => x.Contig != y.Contig && y.Contig.Contains(x.Contig)));

            // Export Data
            SaveData(finalList);
        }

        //private List<ContigModel> SearchContig(ContigModel currentContig, List<ContigModel> route)
        //{
        //    var currentRoute = new List<ContigModel>(route);
        //    if (currentContig.NextContigModels.Count > 1)
        //    {
        //        var toSearch = currentContig.NextContigModels.Keys.OrderByDescending(c => c.NextContigModels.Values.Count).ToList();
        //        foreach (var contigModel in toSearch)
        //        {
        //            if (contigModel.Visited)
        //            {
        //                return null;
        //            }
        //            contigModel.Visited = true;
        //            currentRoute.Add(contigModel);
        //            var a = SearchContig(contigModel, currentRoute);
        //            if (a != null)
        //            {
        //                return a;
        //            }
        //            currentRoute.Remove(contigModel);
        //        }
        //    }
        //    else if (currentContig.NextContigModels.Count == 1)
        //    {
        //        var nextContig = currentContig.NextContigModels.First();
        //        nextContig.Key.Visited = true;
        //        currentRoute.Add(nextContig.Key);
        //        var a = SearchContig(nextContig.Key, currentRoute);
        //        if (a != null)
        //        {
        //            return a;
        //        }
        //        currentRoute.Remove(nextContig.Key);
        //    }
        //    else
        //    {
        //        if (currentContig.VisitedEndNode)
        //        {
        //            return null;
        //        }
        //        currentContig.VisitedEndNode = true;
        //        route.Add(currentContig);
        //    }
        //    return route;
        //}

        private void SaveData(List<ContigModel> contigsToExport)
        {
            List<string> exportStrings = new List<string>();
            foreach (var contigModel in contigsToExport)
            {
                exportStrings.Add($"> Built with k-mer of size {_kMerSize}.");
                exportStrings.Add(contigModel.Contig);
            }
            System.IO.File.WriteAllLines(Directory.GetCurrentDirectory()+$"\\rand.n.n50_{CalculateN50(contigsToExport)}.fa", exportStrings);
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
