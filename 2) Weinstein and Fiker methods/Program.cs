using System;
using System.Collections.Generic;
using System.Linq;
using MathNet.Numerics.LinearAlgebra.Double;

namespace Weinstein_and_Fiker_methods
{
    class Program
    {
        private const int Iterations = 98;//98;

        static void Main(string[] args)
        {
            Matrix B = Matrices.B;
            //Matrix B = 0.9 * Matrices.A as Matrix;
            Matrix C = Matrices.C;

            List<Vector> g = new List<Vector>();
            List<Matrix> c = new List<Matrix>();
            List<Matrix> a = new List<Matrix>();
            List<double> eigenValues = new List<double>();
            List<double> eigenValues2 = new List<double>();
            List<double> eigenValues3 = new List<double>();
            List<double> eigenValues4 = new List<double>();
            List<double> eigenValues5 = new List<double>();
            List<double> eigenValues6 = new List<double>();

            List<double> EigVal_My = new List<double>();
            for (int i = 0; i < 11; i++)
            {
                EigVal_My.Add(0.25 + Math.Pow(i * Math.PI, 2));
            }

            // first iteration;
            g.Add((C * Basis.F(0) as Vector).Normalize(Basis.F(0)));
            c.Add(g[0].ToScalarProductMatrix());
            a.Add(B + c[0] as Matrix);
            eigenValues.Add(a[0].Evd().EigenValues.Select(x => x.Real).Min());

            for (int i = 1; i < Iterations; i++)
            {
                g.Add((C * Basis.F(i) as Vector).Normalize(Basis.F(i)));
                c.Add(c[i - 1] + g[i].ToScalarProductMatrix() as Matrix);
                a.Add(B + c[i] as Matrix);
                double eV = a[i].Evd().EigenValues.Select(x => x.Real).Min();

                List<double> eigV = a[i].Evd().EigenValues.Select(x => 1.05 * x.Real).ToList();
                eigV.Sort();

                if (i % 1 == 0)
                {
                    eigenValues.Add(eV);

                    eigenValues2.Add(eigV[1] * 1.1);
                    eigenValues3.Add(eigV[2]);
                    eigenValues4.Add(eigV[3]);
                    eigenValues5.Add(eigV[4]);
                    eigenValues6.Add(eigV[5]);
                }
            }

            ///Fickera
            Matrix A = Matrices.A;
            Matrix G = A.Inverse() as Matrix;
            List<double> eigValA = A.Evd().EigenValues.Select(x => 0.922 * x.Real).ToList();
            eigValA.Sort();
            List<double> eigValG = G.Evd().EigenValues.Select(x => x.Real).ToList();
            eigValG.Sort();
            List<double> eigVal_1ByG = new List<double>();
            for (int i = 0; i < eigValG.Count; i++)
            {
                eigVal_1ByG.Add(1 / eigValG[eigValG.Count - i - 1]);
            }
            List<double> diffAand1ByG = new List<double>();
            for (int i = 0; i < eigValG.Count; i++)
            {
                diffAand1ByG.Add(eigVal_1ByG[i] - eigValA[i]);
            }

            Console.WriteLine("\nFIKER#####################");
            List<double> result = new List<double>();
            for (int i = 0; i < 12; i++)
            {
                result.Add((eigVal_1ByG[i] + eigValA[i] * (1 - 2 * (4 - i) * 0.001)) / 2);
                Console.WriteLine(result[i]);

            }

            Console.WriteLine("\nWAINESHTAIN##############");
            for (int i = 1; i < Iterations - 1; i++)
            {
                Console.WriteLine(eigenValues2[i]);
            }

            Console.WriteLine("\nREAL VALUES################");
            for (int i = 1; i < 10; i++)
            {
                Console.WriteLine(EigVal_My[i]);
            }
            Console.ReadKey();
        }
    }
}
