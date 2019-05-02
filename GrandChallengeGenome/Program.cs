using System;
using System.Collections.Generic;
using System.Reflection.PortableExecutable;
using GrandChallengeGenome.Models;
using REvernus.Core.Serialization;

namespace GrandChallengeGenome
{
    internal class Program
    {
        private static void Main(string[] args)
        {
            // try and serialize file provided in args
            var data = Serializer.DeserializeData(args[0]);

            var g = new DeBruijn(data);
            g.AssembleGenome(int.Parse(args[1]));
        }
    }
}
