/* Este programa resuelve las ecuaciónes de Euler-Lagrange para el 
 * pendulo elastico por medio del metodo de Runge-Kutta 4
 * El usuario lo debe introducir las condiciones iniciales 
 * y como resultado obtendra un archivo con los datos del resultado */
 
#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <stdlib.h>
#include <string.h>
// Es importante definir el paso y otras variables que dependen de el
#define h 0.0041  		// paso en segundos
#define tt 3700 			// numero de pasos

// Primero definiremos unas constantes importantes:
double g = 9.8; 		// Este es el valor de la aceleración en m/s^2
float kr = 25; 			// Esta es la constante del resorte en N/m
float m = 5; 			// Esta es la masa en kg
float l = 1; 			// Esta es la longitud inicial de resorte en m
double pi = 3.141592; 	// Esta es la constante pi
double M[tt+1][7] ={0};	// Esta es la matriz de variables.
int i;					// Indices generales

// Guardamos las funciones a utilizar:
double dv(double z1, double th1, double w1);
double dw(double z1, double vz1, double th1, double w1);
double RK4(double x, double dx, double y, double dy, double t);

// Aqui empieza el programa:
int main(int u, char *ent[]){
	
	char a[20], b[20], c[20], d[20], f[20], g[20], j[20], n[20], o[20];
	
	// Entrada del programa final03.exe z th kr m l
	// Constantes
	//double e = atof(ent[3]);
	if(ent[3] != 0){
		kr = atof(ent[3]);
		m = atof(ent[4]);
		l = atof(ent[5]);
		strcpy(n,ent[3]);
	}
	else{
		strcpy(n,"def");
	}
	
	// Condiciones iniciales:
	double z, vz, th, w, t;
	t = 0.0;			// tiempo inicial
	z = atof(ent[1]);		// deformacion inicial
	vz = 0.0;			// velocidad linial inicial
	th = (atof(ent[2]))*pi;		// angulo inicial
	w = 0.0;			// velocidad angular inicial
	
	//nombre del archivo de salida:
   strcpy(a,ent[1]);
   strcpy(b,ent[2]);

   strcpy(c,"datos");
   strcpy(d,".txt");
   strcat(a,b);
   strcpy(f,a);
   strcat(f,n);
   strcpy(o,f);
   strcat(c,o);
   strcpy(g,c);
   strcat(g,d);
   strcpy(j,g);
	
	RK4(z,vz,th,w,t);
	
	FILE *datos = fopen(j,"w");
	fprintf(datos, "---- Félix Cabrera, ECFM-USAC, 2019 ---- \n");
	fprintf(datos, "Solución numérica para el péndulo elástico: \n");
	fprintf(datos, "Pendulo de longitud %.2f [m], constante %.2f [N/m] y masa %.2f [kg]. \n",l,kr,m);
	fprintf(datos, " z[m], vz[m/s], theta[rad], w[rad/s], t[s], x[m], y[m], \n");
	for(i=0;i<tt;i++){
		fprintf(datos, "%f, %f, %f, %f, %f, %f, %f, \n",M[i][0],M[i][1],M[i][2],M[i][3],M[i][4],M[i][5],M[i][6]);
	}
	fclose(datos);
	
	printf("Se ha generado un archivo %s con las souciones. \n",j);
	printf("con las condiciones iniciales z = %.2f [m] y  th = %.2f  [rad]\n", z, th);
	printf("Pendulo de longitud %.2f [m], constante %.2f [N/m] y masa %.2f [kg]. \n",l,kr,m);
	printf("utilizando el metodo RK4 con un paso de 0.004 s hasta 15 s. \n");
	return 0;
}
// Definimos las funciones:
double dv(double z1, double th1, double w1){
	double dv = (z1+l)*(w1*w1)+g*cos(th1)-(kr/m)*z1;
	return dv;
}
double dw(double z1, double vz1, double th1, double w1){
	double dw = -(1/(l+z1))*(2*w1*vz1 + g*sin(th1));
	return dw;
}
double RK4(double x, double dx, double y, double dy, double t){
	
	float k1,k2,k3,k4,l1,l2,l3,l4,q1,q2,q3,q4,m1,m2,m3,m4; //Estas son las variables para RK4
	
	M[0][0] = x; // z
	M[0][1] = dx; // dz
	M[0][2] = y; // th
	M[0][3] = dy; // w
	M[0][4] = t;
	
	for(i=0;i<tt;i++){

		//Aqui se aplica el metodo:
		k1 = h*M[i][1];
		q1 = h*M[i][3];
		l1 = h*(dv(M[i][0],M[i][2],M[i][3]));
		m1 = h*(dw(M[i][0],M[i][1],M[i][2],M[i][3]));
		
		k2 = h*(M[i][1]+0.5*l1);
		q2 = h*(M[i][3]+0.5*m1);
		l2 = h*(dv(M[i][0]+0.5*k1,M[i][2]+0.5*q1,M[i][3]+0.5*m1));
		m2 = h*(dw(M[i][0]+0.5*k1,M[i][1]+0.5*l1,M[i][2]+0.5*q1,M[i][3]+0.5*m1));
	
		k3 = h*(M[i][1]+0.5*l2);
		q3 = h*(M[i][3]+0.5*m2);
		l3 = h*(dv(M[i][0]+0.5*k2,M[i][2]+0.5*q2,M[i][3]+0.5*m2));
		m3 = h*(dw(M[i][0]+0.5*k2,M[i][1]+0.5*l2,M[i][2]+0.5*q2,M[i][3]+0.5*m2));
	
		k4 = h*(M[i][1]+l3);
		q4 = h*(M[i][3]+m3);
		l4 = h*(dv(M[i][0]+k3,M[i][2]+q3,M[i][3]+m3));
		m4 = h*(dw(M[i][0]+k3,M[i][1]+l3,M[i][2]+q3,M[i][3]+m3));
		
		M[i+1][0] = (M[i][0] + (0.166666666)*(k1+2*k2+2*k3+k4));
		M[i+1][1] = (M[i][1] + (0.166666666)*(l1+2*l2+2*l3+l4));
		M[i+1][2] = (M[i][2] + (0.166666666)*(q1+2*q2+2*q3+q4));
		M[i+1][3] = (M[i][3] + (0.166666666)*(m1+2*m2+2*m3+m4));
		M[i+1][4] = M[i][4] + h;
		M[i][5] = (l+M[i][0])*sin(M[i][2]);
		M[i][6] = -(l+M[i][0])*cos(M[i][2]);
				
	}
	return 0.0;
}