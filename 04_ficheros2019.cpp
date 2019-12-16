# include <stdio.h>
# include <math.h>
# include <string.h>

struct parametros{
    double alpha, beta, Gamma, omega;
};
typedef struct parametros Parametros;
/** Par�metros del p�ndulo*/

struct estado {
    double theta, w;
};
typedef struct estado Estado;
/** El estado del p�ndulo queda definido
    por los valores del �ngulo y la velocidad angular.
    Suponemos que vienen expresados,respectivamente, en radianes y radianes/segundo  */

struct datosSimulacion{
    double factorDeltaT, factorTCout, factorTStat, factorTSimul;
    Estado estadoInicial;
    Parametros param;
};
typedef struct datosSimulacion DatosSimulacion;
/** Una simulaci�n se define por su estado inicial (�ngulo y velocidad),
    sus par�metros y los factores temporales.*/


/**
ENTRADA: un real, theta (un �ngulo en radianes)
SALIDA: un real: el �ngulo normalizado en el intervalo [-PI, PI):
        Esto es: un �ngulo theta' tal que:
        i)  -PI <= theta' < PI
        ii) theta - theta' es un m�ltiplo entero de 2 PI
*/
double normaliza (double theta);

/**
ENTRADA: cuatro reales: el coeficiente alpha, un �ngulo theta, una velocidad angular w y un valor para la precisi�n epsilon
SALIDA:  dos reales: el �ngulo l�mite y un T caracter�stico (en caso de movimiento peri�dico, el periodo) (por referencia)
         un entero que caracteriza el tipo de soluci�n:
            -1 si el p�ndulo no se voltea
             0 si acaba en la posici�n de equilibrio inestable theta = PI
             1 si tiene energ�a suficiente para voltearse
*/
int analisis (double alpha, double theta, double w, double epsilon,  double & thetaM, double &T);

/**
ENTRADAS: un conjunto de par�metros, el estado del p�ndulo en el instante t, dos reales( el propio t y DeltaT)
SALIDA:  el estado del p�ndulo en el instante t + DeltaT
*/
Estado Runge_Kutta (Parametros param, Estado e, double t, double deltaT);


/** Desarrolla el siguiente subalgoritmo SIN MODIFICAR LA CABECERA.
La funci�n Runge_Kutta hace uso de �l para calcular la evoluci�n del estado del p�ndulo*/
double aceleracion (Parametros param,Estado e, double t){

}


int main(){



    return 0;
}


double normaliza (double theta)
{
    double n =  floor ( (theta+M_PI)/2/M_PI);
    theta -= 2*n*M_PI;

    return theta;
}

int analisis (double alpha, double theta0, double w0, double epsilon,
               double &thetaM, double &T)
{
    int k;
    double E = w0*w0/2/alpha - cos(theta0);

    if (E == 1) {
        double D = tan (normaliza(theta0)/4 + M_PI/4);

        if (D == 0) T = 10/sqrt(alpha); // Damos un valor arbitrario,
                                        //  para evitar dividir por 0 en la siguiente expresi�n
        else
            T = 1/sqrt(alpha) * log ( tan (0.999*M_PI/4 + M_PI/4) / D);

        if (D == 0 ||w0 < 0) thetaM = - M_PI;
        else thetaM = M_PI;
        k = 0;
    }
    else {
            double T0, r;
            if (E < 1) {
                thetaM = acos (-E);
                T0 = 2*M_PI /sqrt(alpha);
                r = sin (thetaM/2);
                k = -1;
            }
            else{
                double A = sqrt( alpha*(1+E)/2 );
                T0 = M_PI/A;
                r = sqrt(alpha)/A;
                thetaM = M_PI;
                k = 1;
            }
            long long N = (long long) ceil (log (epsilon*(1-r*r))/2/log(r));
            double suma = 1, termino=1;
            for (long long n=1; n<= N ; n++){
                termino = termino * pow ( (2*n-1)/2./n *r , 2);
                suma += termino;
            }
            T = T0 * suma;
    }

    return k;
}

Estado Runge_Kutta (Parametros param, Estado e0, double t, double deltaT){
    Estado eAux, eF;
    double divisor[] = {6,3,3,6}, f[]={0,0.5,0.5,1};
    double dTheta, dw;

    eF =  e0;
    for (int i=0; i< 4; i++) {
        eAux.theta = e0.theta + dTheta * f[i];
        eAux.w = e0.w + dw * f[i];

        dTheta = eAux.w * deltaT;
        dw = aceleracion (param, eAux, t + f[i] *deltaT) * deltaT;

        eF.theta += dTheta / divisor[i];
        eF.w += dw / divisor[i];
    }

    return eF;
}
