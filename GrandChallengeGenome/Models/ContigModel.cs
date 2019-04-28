using System;
using System.Collections.Generic;
using System.Text;

namespace GrandChallengeGenome.Models
{
    class ContigModel
    {
        public string Contig { get; set; }

        public HashSet<ContigModel> InContigModels { get; set; }

        public HashSet<ContigModel> OutContigModels { get; set; }
    }
}
