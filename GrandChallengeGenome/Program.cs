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
            var data = Serializer.DeserializeData("C:\\Users\\Connor\\Downloads\\rand.500.1.fq");

            //DoTestData(data); return;

            // loop and build several debrujin graphs
            for (int i = 30; i < 500; i+=20)
            {
                var g = new DeBruijn(data);
                g.BuildDeBruijnGraph(i);
                g.DisposeMemory();
                GC.Collect();
                GC.WaitForPendingFinalizers();
            }
        }

        private static void DoTestData(List<BaseContigModel> data)
        {
            var g = new DeBruijn(data);
            g.BuildDeBruijnGraph(1);
            g.DisposeMemory();
            GC.Collect();
            GC.WaitForPendingFinalizers();
        }
    }
}
