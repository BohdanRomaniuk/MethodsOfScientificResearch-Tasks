using MathNet.Numerics.LinearAlgebra.Double;
using System;

namespace Weinstein_and_Fiker_methods
{
    public static class VectorExtensions
    {
        public static Vector Normalize(this Vector vector, Vector basis)
        {
            //double coef = vector*basis;
            //return vector / Math.Sqrt(coef) as Vector;
            /*
             double coef = vector.AbsoluteMaximum() ;
             return vector / coef as Vector;
 */

            double coef = vector * vector;
            return vector / Math.Sqrt(coef) as Vector;
        }

        public static Vector Trim(this Vector vector)
        {
            Vector result = (Vector)vector.Clone();
            result[0] = 0;
            result[1] = result[2] / 2;
            result[Constants.N - 2] = result[Constants.N - 3] / 2;
            result[Constants.N - 1] = 0;
            return result;
        }

        public static Matrix ToScalarProductMatrix(this Vector vector)
        {
            Matrix matrix = new DenseMatrix(Constants.N, Constants.N);
            for (int i = 0; i < Constants.N; i++)
            {
                for (int j = 0; j < Constants.N; j++)
                {
                    matrix[i, j] = vector[i] * vector[j];
                }
            }
            return matrix;
        }
    }
}
