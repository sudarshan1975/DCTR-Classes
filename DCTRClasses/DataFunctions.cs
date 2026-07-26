namespace DCTRClasses
{
    public class DataFunctions
    {
        public static double[,] ReadArrayFromFile(string fileName)
        {
            List<List<double>> outList = [];

            using (StreamReader sr = new(fileName))
            {
                string? line = sr.ReadLine();

                while (line != null)
                {
                    var lineList = line.Split('\t');

                    List<double> list = [.. lineList.Select(double.Parse)];

                    outList.Add(list);

                    line = sr.ReadLine();
                }
            }

            if (outList.Count == 0) return new double[0, 0];

            int M = outList.Count, N = outList[0].Count;

            var outArray = new double[M, N];

            for (int i = 0; i < M; i++)
            {
                for (int j = 0; j < N; j++)
                {
                    outArray[i, j] = outList[i][j];
                }
            }

            return outArray;
        }
    }
}
