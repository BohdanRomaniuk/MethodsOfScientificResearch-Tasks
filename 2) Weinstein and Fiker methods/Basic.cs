using System;
using MathNet.Numerics.LinearAlgebra.Double;

namespace Weinstein_and_Fiker_methods
{
    public static class Basis
    {
        public static Vector F(int k)
        {
            Vector vector = new DenseVector(Constants.N);
            //vector[k] = 1;

            for (int i = 0; i < Constants.N; i++)
            {
                double x = Constants.Left + i * Constants.H;
                vector[i] = Math.Pow(x, (k));
                //vector[i] = Math.Sin((k + 1) * Math.PI * x);
            }

            return vector;
        }
    }
}
