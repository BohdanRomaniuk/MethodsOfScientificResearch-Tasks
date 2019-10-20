using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace _2__Weinstein_and_Fiker_methods
{
    public class Program
    {
        #region Fields

        public double a { get; set; }
        public double b { get; set; }
        public int number { get; set; }
        public double h { get { return (b - a) / (number - 1); } } // If we have 50 nodes, then we have 49 segments

        public double[,] A { get; set; }
        public double[,] B { get; set; }

        #endregion

        #region Some methods
        /// <summary>
        ///
        /// </summary>
        /// <param name="a">Lower bound</param>
        /// <param name="b">Upper bound</param>
        /// <param name="n">A & B matrix size</param>
        public Program(double a, double b, int n)
        {
            this.A = A;
            this.B = B;
            this.a = a;
            this.b = b;
            this.number = n;
            SetAandB();
        }

        public void SetAandB()
        {
            // T - matrix with only 1/h^2 * [ ... 1 -2 1 ... ]
            double[,] T = new double[number, number];

            #region Only matrix
            T[0, 0] = 1 / Math.Pow(h, 2);
            for (int i = 1; i < number; ++i)
            {
                for (int j = 0; j < number - 1; ++j)
                {
                    if (j == i)
                    {
                        T[i, j] = -2 / Math.Pow(h, 2);
                    }
                    if (j == i - 1 || j == i + 1)
                    {
                        T[i, j] = 1 / Math.Pow(h, 2);
                    }
                }
            }
            T[number - 1, number - 1] = 1 / Math.Pow(h, 2);
            #endregion

            double[,] A = (double[,])T.Clone();
            double[,] B = (double[,])T.Clone();

            for (int i = 0; i < number; ++i)
            {
                double x = Math.Pow(4 - Math.Pow(a + i * h, 2), 2);
                double y = 2 * Math.Pow(4 - Math.Pow(a + i * h, 2), 2);

                for (int j = 0; j < number; ++j)
                {
                    A[i, j] *= -x;
                    B[i, j] *= y;
                }
            }

            this.A = A;
            this.B = B;
        }

        /// <summary>
        /// Return matrix C. C = A - B.
        /// </summary>
        /// <param name="aMatr">A matrix</param>
        /// <param name="bMatr">B matrix</param>
        /// <param name="n">A & B matrix size</param>
        /// <returns></returns>
        public double[,] GetC(double[,] aMatr, double[,] bMatr, int n)
        {
            double[,] cMatr = new double[n, n];
            for (int i = 0; i < n; ++i)
            {
                for (int j = 0; j < n; ++j)
                {
                    cMatr[i, j] = aMatr[i, j] - bMatr[i, j];
                }
            }
            return cMatr;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="k">The power of function</param>
        /// <returns></returns>
        public double[] GetFk(int k)
        {
            double[] fk = new double[number];
            for (int i = 0; i < number; ++i)
            {
                //fk[i] = Math.Pow(a + i * h, k); // x^k
                fk[i] = Math.Sin(k * (a + i * h)); // sin(kx)
            }
            return fk;
        }

        public double[,] GetCn(List<double[]> g, int n)
        {
            double[,] cn = new double[n, n];
            for (int i = 0; i < n; ++i)
            {
                double elem = 0;
                for (int j = 0; j < n; ++j)
                {
                    double b = 0; // TODO: sum(g)/n
                    for (int jj = 0; jj < n; ++jj)
                    {
                        b += g[i][jj];
                    }
                    b /= n;

                    elem += b * g[i][j];
                }
            }
            return cn;
        }

        public double[,] GetAn(double[,] B, double[,] Cn)
        {
            double[,] result = new double[number, number];

            for (int i = 0; i < number; ++i)
            {
                for (int j = 0; j < number; ++j)
                {
                    result[i, j] = B[i, j] + Cn[i, j];
                }
            }

            return result;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="fk">fk</param>
        /// <returns>gk</returns>
        public double[] GetGk(double[,] cMatr, double[] fk)
        {
            // gk = C * fk
            double[] gk = new double[number];
            for (int i = 0; i < number; ++i)
            {
                gk[i] = 0;
                for (int j = 0; j < number; ++j)
                {
                    gk[i] += cMatr[i, j] * fk[j];
                }
            }
            return gk;
        }

        #endregion

        public void Calculate()
        {
            double[,] C = GetC(A, B, number);

            List<double[]> listFk = new List<double[]>();
            for (int k = 0; k < 50; ++k)
            {
                listFk.Add(GetFk(k));
            }

            List<double[]> listGk = new List<double[]>();
            for (int k = 0; k < 50; ++k)
            {
                listGk.Add(GetGk(C, listFk[k]));
            }

            // Find Cn, find An, using lib find ...values of An.
            double[,] Cn = GetCn(listGk, number);
            double[,] An = GetAn(B, Cn);

            double[] d = new double[number];
            double[] e = new double[number];
            double[] result = new double[number];

            for (int i = 0; i < number; ++i)
            {
                d[i] = An[i, i];
                if (i == 0)
                {
                    e[i] = d[i + 1];
                }
                else
                {
                    e[i] = d[i - 1];
                }
            }

            alglib.smatrixtdevd(ref d, e, number, 0, ref An);

            Console.WriteLine(d.Length + " " + e.Length + " " + An.Length);

            for (int i = 0; i < number; i++)
            {
                for (int j = 0; j < number; j++)
                {
                    Console.WriteLine(d[i] + " " + e[i]);
                }
            }
        }
    }
}
