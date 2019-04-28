using System;
using GrandChallengeGenome.Models;
using REvernus.Core.Serialization;

namespace GrandChallengeGenome
{
    internal class Program
    {
        private static void Main(string[] args)
        {
            // try and serialize file provided in args
            var a = Serializer.DeserializeData("C:\\Users\\Connor\\Downloads\\rand.500.1.fq");
            var d = new DeBruijn(a);
            d.BuildDeBruijnGraph(50);
            Console.WriteLine(a);
        }
    }
}
