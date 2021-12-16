using System;

namespace CP3
{
    class Program
    {
        //static int N = 700, M = N;
        static double ax = -10, bx = 10, ay = 0, height =0.5, by = height;
        static double  epsilon = 0.1, Mach = Math.Pow(0.6, 2);

        static double[,] Relaxation(double[] x, double[] y, int N, int M, double hx, double hy)
        {
            int counter = 0;
            double gamma = Math.Pow(hy / hx, 2);
            double[,] previous = new double[N, M], psi = new double[N, M];
            bool flag = true;
            //double omega = 0.1;
            
             for (int j = 1; j < M ; j++)
             {
                psi[0, j] = 0;
                psi[N - 1, j] = 0;
             }

             for (int i = 0; i < N-1; i++)
             {
                 psi[i, 0] = 0;
                 psi[i, M - 1] = 0;
                 if (x[i] >= -1 && x[i] <= 1)
                    psi[i, 0] =  -epsilon* (1 - Math.Pow(x[i], 2));
                 
             }
             
             for (int i = 0; i < N; i++)
                for (int j = 0; j < M; j++)
                  previous[i, j] = psi[i, j];
             



             while (flag)
             {
                 int k = 0;
                 for (int i = 1; i < N - 1; i++)
                     for (int j = 1; j < M - 1; j++)
                     {
                        //r[i, j] = (1 - Mach) * gamma * (previous[i + 1, j] - 2 * previous[i, j] + previous[i - 1, j]) + previous[i, j + 1] - 2 * previous[i, j] + previous[i, j - 1];
                        psi[i, j] = 1/(2 * (gamma * (1 - Mach) + 1)) * (gamma * (1 - Mach) * (previous[i - 1, j] + previous[i + 1, j]) + previous[i, j - 1] + previous[i, j + 1]);
                        //psi[i, j] = omega*psi[i, j] +  (1-omega)*previous[i,j];
                        //psi[i, j] = previous[i, j] + omega * r[i, j];
                     }

                 for (int i = 0; i < N; i++)
                     for (int j = 0; j < M; j++)
                         if (Math.Abs(psi[i, j] - previous[i, j]) < (10e-5))
                             k++;

                 if (k == M * N)
                     flag = false;


                 for (int i = 0; i < N; i++)
                     for (int j = 0; j < M; j++)
                         previous[i, j] = psi[i, j];
                counter++;

                Console.WriteLine("{0};{1}", counter, 2 * (1 / (1 - Mach) * (psi[N / 2, 2] - psi[N / 2, 0]) / (2 * hy)));
            }
            //Console.WriteLine("{0}; {1}", N, 2 * (1 / (1 - Mach) * (psi[N / 2, 2] - psi[N / 2, 0]) / (2 * hy)));
            //Console.WriteLine("{0}; {1}", N, 2 * (1 / (1 - Mach) * (psi[N/2, 2] - psi[N/2, 0]) / (2 * hy)));
             //for (int i = 0; i < N; i++)
             // Console.WriteLine("{0}; {1}", x[i], 2 * (1 / (1 - Mach) * (psi[i, 2] - psi[i, 0]) / (2 * hy)));
            return psi;
        }
       

        static double[,] SIP(double[] x, double[] y, int N, int M, double hx, double hy) {
           double gamma = Math.Pow(hy / hx, 2);
            double[,] previous = new double[N, M], psi = new double[N, M];
            bool flag = true;
            int counter = 0, k;
           
            double[] A = new double[M], B = new double[M], C = new double[M], F = new double[M], alpha = new double[M], beta = new double[M];
            for (int i = 0; i < N; i++)
            {
                for (int j = 0; j < M; j++)
                {
                    psi[i, j] = 0;
                    previous[i, j] = 0;
                }
            }




            while (flag)
            {
                k = 0;
                for (int i = 1; i < N - 1; i++)
                {

                    A[0] = 0;
                    B[0] = 1;
                    F[0] = 0;
                    if (x[i] >= -1 && x[i] <= 1)
                        F[0] = -epsilon * (1 - Math.Pow(x[i], 2));

                    B[M - 1] = 1;
                    C[M - 1] = 0;
                    F[M - 1] = 0;

                    alpha[1] = -C[0] / B[0];
                    beta[1] = F[0] / B[0];

                    for (int j = 1; j < M - 1; j++)
                    {
                        A[j] = 1;
                        B[j] = -2 * (gamma * (1 - Mach) + 1);
                        C[j] = 1;
                        F[j] = -gamma * (1 - Mach) * (previous[i + 1, j] + psi[i - 1, j]);

                        alpha[j + 1] = -C[j] / (A[j] * alpha[j] + B[j]);
                        beta[j + 1] = (F[j] - A[j] * beta[j]) / (A[j] * alpha[j] + B[j]);

                    }

                    psi[i, M - 1] = (F[M - 1] - A[M - 1] * beta[M - 1]) / (B[M - 1] + A[M - 1] * alpha[M - 1]);

                    for (int j = M - 2; j >= 0; j--)
                    {
                        psi[i, j] = alpha[j + 1] * psi[i, j + 1] + beta[j + 1];

                     
                    }

                }
                for (int i = 0; i < N; i++)
                    for (int j = 0; j < M; j++)
                        if (Math.Abs(psi[i, j] - previous[i, j]) < (10e-6))
                            k++;
                    
                if (k == M * N)
                    flag = false;


                for (int i = 0; i < N; i++)
                    for (int j = 0; j < M; j++)
                        previous[i, j] = psi[i, j];

                
                counter++;

                //Console.WriteLine("{0};{1}", counter, 2 * (1 / (1 - Mach) * (psi[N/2, 2] - psi[N/2, 0]) / (2 * hy)));
            }

            Console.WriteLine("{0}; {1}", hx, 2 * (1 / (1 - Mach) * (psi[N/2, 2] - psi[N/2, 0]) / (2 * hy)));
           // for(int i=0;i<N;i++)
            //Console.WriteLine("{0}; {1}", x[i], 2 * (1 / (1 - Mach) * (psi[i, 2] - psi[i, 0]) / (2 * hy)));

            return psi;
        }
        static void Main(string[] args)
        {

            int N, M;
           

            for (N = 1000; N <=1000; N += 10)
            {
                double hx = (bx - ax) / Convert.ToDouble(N), hy = hx;
                M = Convert.ToInt32((by - ay) / hx);
                double[] x = new double[N], y = new double[M];
               
                x[0] = ax;
                for (int i = 0; i < N - 1; i++)
                    x[i + 1] = x[i] + hx;

                y[0] = ay;
                for (int j = 0; j < M - 1; j++)
                    y[j + 1] = y[j] + hy;
                
                //Relaxation(x, y, N, M, hx, hy);
                SIP(x, y, N, M, hx, hy);
            }
            //

        }

    }
}
