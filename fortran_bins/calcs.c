

//----------------------------------------------------------------------
//----------------------------------------------------------------------
double Slow_Calc_Semblance(strSignal *data,float T, double s, strToolSTC S)
{
    double numerador,denominador;
    double a,p,t,temp,step;
    int i,r,reso;

    // Observação.
    step = 0.1; // resolução no tempo : colocar como parâmetro
    reso = (int) ( (double) S.W / step ) ;

    //--------------------------------------------------------------------
    //--------------------------------------------------------------------
    numerador = 0;
    for ( i = 0; i < reso; i ++ ) // Soma em T
    {
        t = T + i * step;
        a = 0;
        for ( r = 0; r < S.nD; r++) // Soma em M
        {
            p = t + (s * r * S.D)/S.TA;
            if ( p < S.nT && p >= 0) // Numerador
            {

                temp = LI(data,r,p);
                a = a + temp;
            }
        }
        numerador = numerador + a * a;
    }
    //--------------------------------------------------------------------
    //--------------------------------------------------------------------
    denominador = 0;
    for ( i = 0; i < reso; i ++ )
    {
        a = 0;
        t = T + i * step;
        for ( r = 0; r < S.nD; r++)
        {
            p = t + (s * r * S.D)/S.TA;
            if ( p < S.nT && p >= 0)
            {
                temp = LI(data,r,p);
                a = a + temp*temp;
            }
        }
        denominador = denominador + a;
    }
    if ( denominador == 0 )
        return 0;
    else
    {
        temp = numerador / denominador;
        return (temp);
    }
}