#include <stdio.h>
#include <math.h>
#include <omp.h>
#include <R.h>
#include <Rdefines.h>

double lfact(long int *a)
{
    int b = 2;
    double result = 0;
    while((*a)>=b)
    {
        result += log(b);
        b++;
    }

    return result;
}


long int min4(long int *a,long int *b,long int *c,long int *d)
{
    long int r = *a;
    if(*b<r)
    {
        r = *b;
    }
    if(*c<r)
    {
        r = *c;
    }
    if(*d<r)
    {
        r = *d;
    }
    return r;
}

void fisher(long int *a, long int *b, long int *c, long int *d, double *la, double *r)
{
    long int t1 = *a+*b;
    long int *p1 = &t1;
    long int t2 = *c+*d;
    long int *p2 = &t2;
    long int t3 = *c+*a;
    long int *p3 = &t3;
    long int t4 = *b+*d;
    long int *p4 = &t4;
    long int n = t1+t2;

    double currentD = fabs((double)(*a)*(*d)-(*b)*(*c));
    //double currentlnP = lfact(p1)+lfact(p2)+lfact(p3)+lfact(p4)-lfact(pn)-lfact(a)-lfact(b)-lfact(c)-lfact(d);
    double currentlnP = *(la+t1)+*(la+t2)+*(la+t3)+*(la+t4)-*(la+n)-*(la+(*a))-*(la+(*b))-*(la+(*c))-*(la+(*d));

    long int condition = min4(p1,p2,p3,p4);
    //printf("%ld %ld %ld %ld\n",*p1,*p2,*p3,*p4);
    //printf("%ld\n",condition);
    long int a_tmp, b_tmp, c_tmp, d_tmp, i;  
    double D,lnP,P;

    printf("%ld %ld %ld %ld     %lf %lf\n",*a,*b,*c,*d,currentD,currentlnP);

    if(condition==*p1||condition==*p3)
    {
        for(i=0;i<=condition;i++)
        {
            a_tmp = i;
            b_tmp = (*p1)-i;
            c_tmp = (*p3)-i;
            d_tmp = (*p4)-b_tmp;

            D = fabs((double)(i*d_tmp-b_tmp*c_tmp));
            if(D>=currentD)
            {
                //lnP = lfact(p1)+lfact(p2)+lfact(p3)+lfact(p4)-lfact(pn)-lfact(a_tmp_p)-lfact(b_tmp_p)-lfact(c_tmp_p)-lfact(d_tmp_p);
                lnP = *(la+t1)+*(la+t2)+*(la+t3)+*(la+t4)-*(la+n)-*(la+(a_tmp))-*(la+(b_tmp))-*(la+(c_tmp))-*(la+(d_tmp));
                printf("1--%ld %ld %ld %ld     %lf %lf\n",a_tmp,b_tmp,c_tmp,d_tmp,D,lnP);

                if(fabs(lnP-currentlnP)<1e-8||lnP<currentlnP)
                {
                    P = exp(lnP);
                    *r += P;
                    printf("%0.32lf\n",P);
                }
            }
        }
    }
    else
    {
        for(i=0;i<=condition;i++)
        {
            d_tmp = i;
            b_tmp = (*p4)-i;
            c_tmp = (*p2)-i;
            a_tmp = (*p1)-b_tmp;

            D = fabs((double)(i*a_tmp-b_tmp*c_tmp));
            printf("2--%ld %ld %ld %ld     %lf ",a_tmp,b_tmp,c_tmp,d_tmp,D);
            if(D>=currentD)
            {
                //lnP = lfact(p1)+lfact(p2)+lfact(p3)+lfact(p4)-lfact(pn)-lfact(a_tmp_p)-lfact(b_tmp_p)-lfact(c_tmp_p)-lfact(d_tmp_p);
                lnP = *(la+t1)+*(la+t2)+*(la+t3)+*(la+t4)-*(la+n)-*(la+a_tmp)-*(la+b_tmp)-*(la+c_tmp)-*(la+d_tmp);
                printf("%0.32lf\n",lnP);
                
                if(fabs(lnP-currentlnP)<1e-8 || lnP<currentlnP)
                {
                    P = exp(lnP);
                    *r += P;
                }
            }
            else{printf("\n");}
        }
    }
}


void FisherExactTest(long int *total, long int *length, double *la ,double *r)
{


    long int i;


    printf("%ld\n", *length);
    #pragma omp parallel for
    for(i=0;i<*length;i+=4)
    {
        fisher(total+i, total+i+1, total+i+2, total+3+i, la, r+(i/4));
    }

}

int main()
{
    double c1[2]={0,0};
    /*long int t1 = 43;
    long int t2 = 50;
    long int t3 = 7;
    long int t4 = 0;
    long int *p1 = &t1;
    long int *p2 = &t2;
    long int *p3 = &t3;
    long int *p4 = &t4;*/

    long int length = 8;
    long int *l = &length;

    long int a[8] = {43,50,7,0,4,15,18,6};
    long int *arr = a;

    double lfact_array[201];
    
    long int *p;
    double *la = lfact_array;
    long int i;

    for(i=0;i<201;i++)
    {
        p = &i;
        lfact_array[i] = lfact(p);
    }

    double *c = &c1[0];
    FisherExactTest(arr, l, la, c);
    printf("%0.32lf",c[0]);
    return 0;
}