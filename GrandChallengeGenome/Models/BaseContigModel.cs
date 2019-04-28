using System;
using System.Collections.Generic;
using System.Text;

namespace GrandChallengeGenome.Models
{
    /// <summary>
    /// Base model for Contigs read in from the .fq file.
    /// We keep this and convert to a normal one later in case we want to add/keep data for later.
    /// </summary>
    public class BaseContigModel
    {
        public string Contig { get; set; }
    }
}
