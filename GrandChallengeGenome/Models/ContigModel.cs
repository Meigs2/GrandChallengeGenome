using System;
using System.Collections.Generic;
using System.Text;

namespace GrandChallengeGenome.Models
{
    class ContigModel
    {
        public string Contig { get; set; }

        public HashSet<ContigModel> PreviousContigModels { get; set; } = new HashSet<ContigModel>();

        public HashSet<ContigModel> NextContigModels { get; set; } = new HashSet<ContigModel>();

        public bool MarkedForDeletion { get; set; } = false;
    }
}
