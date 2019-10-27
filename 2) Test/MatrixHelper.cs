using System;

namespace _2__Test
{
    public static class MatrixHelper
    {
        public static void Print(this double[,] matrix)
        {
            for (int i = 0; i < matrix.GetLength(0); ++i)
            {
                for (int j = 0; j < matrix.GetLength(1); ++j)
                {
                    Console.Write($"{matrix[i, j]}\t");
                }
                Console.WriteLine();
            }
        }

        public static void Print(this double[] vector)
        {
            for (int i = 0; i < vector.Length; ++i)
            {
                Console.Write($"{vector[i]}\t");
            }
            Console.WriteLine();
        }
    }
}
