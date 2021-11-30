using System;

namespace CP3
{
    class Program
    {
        static int N = 1000, M = N;
        static double ax = -10, bx = 10, ay = 0, height = 10, by = height;
        static double hx = (bx - ax) / Convert.ToDouble(N), hy = (by - ay) / Convert.ToDouble(M), epsilon = 0.1, gamma = Math.Pow(hy / hx, 2), Mach = Math.Pow(0.6, 2);

        static double[,] Relaxation(double[] x, double[] y)
        {
            double[,] previous = new double[N, M], psi = new double[N, M], r = new double[N,M];
            bool flag = true;
            int counter = 0;
            double omega = 0.1;
            
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

                 for (int i = 1; i < N - 1; i++)
                     for (int j = 1; j < M - 1; j++)
                         if (Math.Abs(psi[i, j] - previous[i, j]) < (10e-6))
                             k++;

                 if (k == (M - 2) * (N - 2))
                     flag = false;


                 for (int i = 1; i < N-1; i++)
                     for (int j = 1; j < M-1; j++)
                         previous[i, j] = psi[i, j];


                 counter++;
                Console.WriteLine("{0};{1}", counter, 2 * (1 / (1 - Mach) * (psi[N / 2, 2] - psi[N / 2, 0]) / (2 * hy)));
            }

          
           

            //for (int i = 0; i < N; i++)
                //if (x[i] >= -1.5 && x[i] <= 1.5)
                    //Console.WriteLine("{0}; {1}", x[i], 2 * (1 / (1 - Mach) * (psi[i, 2] - psi[i, 0]) / (2 * hy)));
              

            

            return psi;
        }
        /*
        static double[,] SIP(double[] x, double[] y)
        {
            double[,] psi = new double[N, M], previous= new double[N, M], B = new double[N, M], C = new double[N, M], D = new double[N, M], E = new double[N, M], F = new double[N, M], G = new double[N, M], H = new double[N, M], V = new double[N, M], R = new double[N, M],
                q = new double[N,M], b = new double[N, M], c = new double[N, M], d = new double[N, M], e = new double[N, M], f = new double[N, M], delta= new double[N,M];
            double alpha = 0.9;

            /* for (int j = 1; j < M; j++)
              {
                  psi[0, j] = 0;
                  psi[1, j] = 0;
                  psi[N - 1, j] = 0;
              }

              for (int i = 0; i < N; i++)
              {
                  psi[i, 0] = 0;
                  psi[i, M - 1] = -height;
                  if (x[i] >= -1 && x[i] <= 1)
                  {
                      psi[i, 0] = -epsilon * (1 - Math.Pow(x[i], 2));
                  }
              }
            for (int i = 0; i < N; i++)
                for (int j = 1; j < M; j++)
                {
                    psi[i, j] = 0;
                    previous[i, j] = 0;
                }

            for ( int i=1;i<N-1;i++)
                for(int j = 1; j < M-1; j++)
                {
                   
                    B[i, j] = 1;
                    D[i, j] = gamma * (1 - Mach);
                    E[i, j] = -2 * (gamma * (1 - Mach) +1);
                    F[i, j] = gamma * (1 - Mach);
                    H[i, j] = 1;
                    q[i, j] = 0;
                }


            E[0, 0] = 1;
            F[0, 0] = 0;
            H[0, 0] = 0;
            E[N - 1, M - 1] = 1;
            q[N - 1, M - 1] = -height;
            B[N - 1, M - 1] = 0;
            D[N - 1, M - 1] = 0;
        

            d[0, 0] = E[0, 0];
            G[0, 0] = 0;

            E[0, M - 1] = 1;




            for (int i = 0; i < N; i++)
            {

                E[i, 0] = 1;
                D[i, 0] = 0;
                F[i, 0] = 0;
                H[i, 0] = 0;
                q[i, 0] = 0;


                if (x[i] >= 1 && x[i] <= 1)
                    q[i, 0] = -epsilon * (1 - Math.Pow(x[i], 2));

                E[i, M - 1] = 1;
                q[i, M - 1] = -height;
                B[i, M - 1] = 0;
                D[i, M - 1] = 0;
                F[i, M - 1] = 0;

                b[i, 0] = B[i, 0];
               
                c[i, 0] = D[i, 0];
                if (i >= 1)
                {
                    d[i, 0] = E[i, 0] + alpha * c[i, 0] * f[i - 1, 0] - c[i, 0] * e[i - 1, 0];
                    
                    G[i, 0] = c[i, 0] * f[i - 1, 0];
                }
                e[i, 0] = F[i, 0] / d[i, 0];
               
                f[i, 0] = H[i, 0] / d[i, 0];
                C[i, 0] = 0;
               
            }
            for(int j = 1; j < M-1; j++)
            {
                E[0, j] = 1;
                B[0, j] = 0;
                F[0, j] = 0;
                H[0, j] = 0;
                E[N - 1, j] = 1;
                B[N - 1, j] = 0;
                D[N - 1, j] = 0;
                H[N - 1, j] = 0;

                b[0, j] = B[0, j] / (1 + alpha * e[0, j - 1]);
                c[0, j] = D[0, j];
                d[0, j] = E[0, j] + alpha * (b[0, j] * e[0, j - 1]) - b[0, j] * f[0, j - 1] ;
                e[0, j] = (F[0, j] - alpha * b[0, j] * e[0, j - 1]) / d[0, j];
                f[0, j] = H[0, j] / d[0, j];
               
                C[0, j] = b[0, j] * e[0, j - 1];
                G[0, j] = 0;
            }

            for (int i = 1; i < N - 1; i++)
                for (int j = 1; j < M - 1; j++)
                {
                    b[i, j] = B[i, j] / (1 + alpha * e[i, j - 1]);
                   
                    c[i, j] = D[i, j] / (1 + alpha * f[i - 1, j]);
                    d[i, j] = E[i, j] + alpha * (b[i, j] * e[i, j - 1] + c[i, j] * f[i - 1, j]) - b[i, j] * f[i, j - 1] - c[i, j] * e[i - 1, j];
                    e[i, j] = (F[i, j] - alpha * b[i, j] * e[i, j - 1]) / d[i, j];
                    f[i, j] = (H[i, j] - alpha * c[i, j] * f[i - 1, j]) / d[i, j];
                    
                    C[i, j] = b[i, j] * e[i, j - 1];
                    G[i, j] = c[i, j] * f[i - 1, j];
                }
            for (int i = 0; i < N; i++)
                for (int j = 0; j < M; j++)
                {
                   
                    delta[i, j] = 0;
                    V[i, j] = 0;
                }
           

          
            bool flag = true;
            int k = 0, counter = 0;
            while (flag)
            {
                k = 0;
                for (int i = 1; i < N - 1; i++)
                    for (int j = 1; j < M - 1; j++)
                    {
                        R[i, j] = q[i, j] - (B[i, j] * previous[i, j - 1] + D[i, j] * previous[i - 1, j] + E[i, j] * previous[i, j] + F[i, j] * previous[i + 1, j] + H[i, j] * previous[i, j + 1]);
                        Console.WriteLine(R[i, j]);
                    }
                for (int i = 1; i < N; i++)
                    for (int j = 1; j < M; j++)
                    {
                        V[i, j] = 1 / d[i, j] * (R[i, j] - b[i, j] * V[i, j - 1] - c[i, j] * V[i - 1, j]);
                      
                    }


                for (int i = 1; i < N - 1; i++)
                    for (int j = 1; j < M - 1; j++)
                    {
                       delta[i+1, j] = 1 / e[i, j] * (V[i, j] - delta[i, j] - f[i, j] * delta[i, j+1]);
                      

                    }


                for (int i = 0; i < N; i++)
                    for (int j = 0; j < M; j++)
                    {
                        psi[i, j] = previous[i, j] + delta[i, j];
                        Console.WriteLine(delta[i, j]);

                    }
                for (int i = 1; i < N - 1; i++)
                    for (int j = 1; j < M - 1; j++)
                        if (Math.Abs(psi[i, j] - previous[i, j]) < (10e-3))
                            k++;

                if (k == (M - 2) * (N - 2) )
                    flag = false;

               
                counter++;
                Console.WriteLine(counter);
                for (int i = 1; i < N - 1; i++)
                    for (int j = 1; j < M - 1; j++)
                        previous[i, j] = psi[i, j];
              
            }

            for (int i = 1; i < N - 1; i++)
            {
                //Console.WriteLine();
                for (int j = 1; j < M - 1; j++)
                {
                  //  if (x[i] >= -1 && x[i] <= 1)
                       // Console.WriteLine(psi[i, j]);
                }
            }
                    return psi;
        }
*/

        static double[,] SIP(double[] x, double[] y) {
            double[,] previous = new double[N, M], psi = new double[N, M];
            bool flag = true;
            int counter = 0;
            double k = 0;
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
                    {


                        if (Math.Abs(psi[i, j] - previous[i, j]) < (10e-5))
                            k++;
                    }
                if (k == M * N)
                    flag = false;


                for (int i = 0; i < N; i++)
                    for (int j = 0; j < M; j++)
                        previous[i, j] = psi[i, j];

                
                counter++;

                Console.WriteLine("{0};{1}", counter, 2 * (1 / (1 - Mach) * (psi[N/2, 2] - psi[N/2, 0]) / (2 * hy)));
            }
           // double[] CP = new double[N];
          //  double max = 0;

            //for (int i = 0; i < N; i++)
            {
              //  CP[i] = 2 * (1 / (1 - Mach) * (psi[i, 2] - psi[i, 0]) / (2 * hy));
               // Console.WriteLine(CP[i]);
            }
          //  max = CP[0];
            //for (int i = 1; i < N; i++)
                //if (CP[i] > max)
                  //  max = CP[i];
           
                //if (x[i] >= -1.5 && x[i] <= 1.5)
                //Console.WriteLine("{0}; {1}", x[i], 2 * (1 / (1 - Mach) * (psi[i, 2] - psi[i, 0]) / (2 * hy)));




                return psi;
        }
        static void Main(string[] args)
        {
            double[] x = new double[N], y = new double[M];
           
            x[0] = ax;
            for (int i = 0; i < N - 1; i++)
            {
                x[i + 1] = x[i] + hx;
            }
            y[0] = ay;
            for (int j = 0; j < M - 1; j++)
            {
                y[j + 1] = y[j] + hy;
            }

            Relaxation(x, y);
            //SIP(x, y);
           
            

           

        }

    }
}
