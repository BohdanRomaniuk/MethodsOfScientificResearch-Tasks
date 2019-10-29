using System;

namespace WeinsteinMethod
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

        public static double[,] Multiply(this double[,] first, double[,] second)
        {
            int rows = first.GetLength(0);
            int cols = first.GetLength(1);
            var result = new double[rows, cols];
            for (int i = 0; i < rows; ++i)
            {
                for (int j = 0; j < cols; ++j)
                {
                    for (int k = 0; k < cols; ++k)
                    {
                        result[i, j] += first[i, k] * second[k, j];
                    }
                }
            }
            return result;
        }
    }
}
