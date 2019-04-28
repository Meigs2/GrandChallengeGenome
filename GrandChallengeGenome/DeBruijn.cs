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
        private readonly List<BaseContigModel> _baseContigModels;
        private readonly List<ContigModel> _contigModels;

        /// <summary>
        /// This dictionary stores the kMer sized substring and returns the related contig
        /// </summary>
        private Dictionary<string, ContigModel> _contigDictionary = new Dictionary<string, ContigModel>();

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
                model.Contig = model.Contig;
                _contigModels.Add(model);
            }
        }

        /// <summary>
        /// Builds a DeBruijn Graph based on the contigs given in the constructor spaced by a specified k-mer.
        /// </summary>
        /// <param name="kMer"></param>
        public void BuildDeBruijnGraph(int kMer)
        {
            Console.WriteLine($"Building Graph with k-mer of size {kMer}.");

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
                    var kMerContig = _contigDictionary[kMerContigSubstring];

                    // if we're looking at the first kmer substring of larger contig, set previous to current and go to the next substring.
                    if (previousKMerContig == null)
                    {
                        previousKMerContig = kMerContig;
                        continue;
                    }

                    // add contig to dict if not already included, and re-set current contig to make sure current contig is the one in the dictionary
                    if (kMerContig == null)
                    {
                        _contigDictionary.Add(kMerContigSubstring, new ContigModel(){Contig = kMerContigSubstring});
                        kMerContig = _contigDictionary[kMerContigSubstring];
                    }



                }
            }
        }
    }
}
