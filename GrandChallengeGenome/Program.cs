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
            var a = Serializer.DeserializeData<BaseContigModel>("C:\\Users\\Connor\\Downloads\\rand.500.1.fq");
            Console.WriteLine(a);
        }
    }
}
