/* Convert a rank-4 tensor in Absolute notation to Observer notation */
void ElasticityT4AbsoluteToObserver
(M3 observer, T4 absolute_notation, T4 observer_notation)
{
    int i,j,k,l, ip,jp,kp,lp, a,b;
    for (i=0; i<3; i++)
        for (j=i; j<3; j++)
        {
            a = VoigtIndex(i,j);
            for (k=0; k<3; k++)
                for (l=k; l<3; l++)
                {
                    b = VoigtIndex(k,l);
                    observer_notation[a][b] = 0;
                    /* pgcc bug here if "-funroll-loops" flag is turned on */
                    for (ip=0; ip<3; ip++)
                        for (jp=0; jp<3; jp++)
                            for (kp=0; kp<3; kp++)
                                for (lp=0; lp<3; lp++)
                                    observer_notation[a][b] +=
                                        T4Element
                                        (absolute_notation,ip,jp,kp,lp) *
                                        observer[i][ip] *
                                        observer[j][jp] *
                                        observer[k][kp] *
                                        observer[l][lp];
                    observer_notation[b][a] = observer_notation[a][b];
                }
        }
    return;
} /* end ElasticityT4AbsoluteToObserver() */

#ifdef _ElasticityT4AbsoluteToObserver_TEST
/* fcc example in Hirth & Lothe, pp. 434 */
int main (int argc, char *argv[])
{
    double T11=10.21,T12=3.34,T44=5.78;
    double t11,t12,t13,t16,t22,t44,t55,H; /* trigonal */
    Elasticity e = UnknownElasticity;
    H = 2*T44 + T12 - T11;
    t11 = T11 + H / 2;
    t12 = T12 - H / 3;
    t13 = T12 - H / 6;
    t16 = H * sqrt(2.) / 6;
    t22 = T11 + H / 3 * 2;
    t44 = T44 - H / 3;
    t55 = T44 - H / 6;
    /* e.O[0][0] = 1; */
    /* e.O[0][1] = -2; */
    /* e.O[0][2] = 1; */
    /* e.O[1][0] = 1; */
    /* e.O[1][1] = 1; */
    /* e.O[1][2] = 1; */
    e.O[0][0] = 1;
    e.O[0][1] = 0;
    e.O[0][2] = 0;
    e.O[1][0] = 0;
    e.O[1][1] = 1;
    e.O[1][2] = 0;
    ElasticityObserverComplete(e.O);
    S3PR("observer = %M\n ", e.O);
    ElasticityCubicT4Assign(T11,T12,T44,e.C);
    S6PR("C = %M\n ",e.C);
    ElasticityT4AbsoluteToObserver(e.O,e.C,e.c);
    S6PR("c = %M\n ",e.c);
    printf ("t11 = %f\n", t11);
    printf ("t12 = %f\n", t12);
    printf ("t13 = %f\n", t13);
    printf ("t16 = %f\n", t16);
    printf ("t22 = %f\n", t22);
    printf ("t44 = %f\n", t44);
    printf ("t55 = %f\n", t55);
    return (0);
}
#endif /* _ElasticityT4AbsoluteToObserver_TEST */
