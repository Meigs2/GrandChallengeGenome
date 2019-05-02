using System;
using System.Collections.Generic;
using System.Text;

namespace GrandChallengeGenome.Models
{
    public class ContigModel
    {
        public string Contig { get; set; } = "";

        public Dictionary<ContigModel, int> PreviousContigModels { get; set; } = new Dictionary<ContigModel, int>();

        public Dictionary<ContigModel, int> NextContigModels { get; set; } = new Dictionary<ContigModel, int>();



        public bool Visited { get; set; } = false;

        public bool MarkedForDeletion { get; set; } = false;
    }
}
