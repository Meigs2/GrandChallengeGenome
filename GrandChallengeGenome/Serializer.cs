using System;
using System.Collections.Generic;
using System.IO;
using System.Runtime.Serialization;
using System.Runtime.Serialization.Formatters.Binary;
using GrandChallengeGenome.Models;

namespace REvernus.Core.Serialization
{
    public static class Serializer
    {
        /// <summary>
        /// Tries to serialize an object to a file provided. If the file path does not exist, it will create it.
        /// </summary>
        /// <typeparam name="T"></typeparam>
        /// <param name="data"></param>
        /// <param name="path"></param>
        public static void SerializeData<T>(T data, string path)
        {
            try
            {
                // Check if the parent directory of the file exist, if not, create it.
                if (!Directory.Exists(Path.GetDirectoryName(path)))
                {
                    Directory.CreateDirectory(Path.GetDirectoryName(path));
                }

                var fileStream = new FileStream(path, FileMode.Create);

                if (data != null)
                {
                    BinaryFormatter formatter = new BinaryFormatter();
                    formatter.Serialize(fileStream, data);
                }
            }
            catch (Exception e)
            {
                Console.WriteLine(e);
                throw;
            }
        }

        /// <summary>
        /// Tries to deserialize data into the type given, in this case will be out BaseContigModel
        /// </summary>
        /// <typeparam name="T"></typeparam>
        /// <param name="path"></param>
        /// <returns></returns>
        public static List<BaseContigModel> DeserializeData(string path)
        {
            try
            {
                var fileStream = new FileStream(path, FileMode.Open);
                StreamReader sr = new StreamReader(fileStream);
                var contigModelList = new List<BaseContigModel>();
                while (!sr.EndOfStream)
                {
                    var model = new BaseContigModel();
                    sr.ReadLine();
                    model.Contig = sr.ReadLine();
                    sr.ReadLine();
                    sr.ReadLine();
                    contigModelList.Add(model);
                }
                fileStream.Close();
                return contigModelList;
            }
            catch (Exception)
            {
                return null;
            }
        }
    }
}
