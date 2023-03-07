#include<iostream>
#include<bits/stdc++.h>
using namespace std;

bool isDDM(double m[][], int n)
{
    for (int i = 0; i < n; i++)
   {       
        double sum = 0;
        for (int j = 0; j < n; j++)            
            sum += abs(m[i][j]);       
 
        sum -= abs(m[i][i]);
 
        if (abs(m[i][i]) < sum)
            return false;
        
    }
 
    return true;
}

int main(){
    int m, iter;
    cout<<"Enter number of equations: ";
    cin>>m;

    // matrix A
    cout<<"Enter A:\n";
    double A[m][m];
    for(int i = 0; i < m ; i++){
        for(int j = 0; j < m; j++){
            cout<<"A["<<i+1<<", "<<j+1<<"] =";
            cin>>A[i][j];
        }
    }

    // matrix b
    cout<<"\nEnter B:\n";
    double b[m];
    for(int i = 0; i < m; i++){
        cout<<"B["<<i+1<<"] :";
        cin>>b[i];
    }

    // solution array and initial assumption
    double x[m], x_prev[m], x_temp[m];
    for(int i = 0; i < m; i++){
        x[i] = 0;
        x_prev[i] = 0;
        x_temp[i] = 0;
    }

    if(isDDM(A,m) == false){
        cout<<"Matrix is not diagonally dominant"<<endl;
        return 1;
    }

    //number of iterations
    cout<<"\nEnter number of iterations :";
    cin>>iter;

    bool equal = false;
    for(int i = 0; i < iter; i++){
        
        // gauss seidel
        for(int p = 0; p < m; p++){
            x[p] = (b[p]/A[p][p]);
            for(int q = 0; q < m; q++){
                if(p == q)
                    continue;
                x[p] = x[p] - ((A[p][q]/A[p][p]) * x_temp[q]);
                x_temp[p] = x[p];
            }
        }

        //print solution
        cout<<"Iteration "<<(i+1)<<" :\n   ";
        for(int p = 0; p < m; p++){
            printf("x%d = %f    ", p + 1, x[p]);
        }
        cout<<endl;

        //kill on equal
        if(i != 0){
            equal = true;
            for(int p = 0; p < m; p++){
                if(round(x[p]*1000000)/1000000 != round(x_prev[p]*1000000)/1000000){
                    equal = false;
                    break;
                }
            }
            if(equal == true)
                break;
        }

        for(int p = 0; p < m; p++){
            x_prev[p] = x[p];
        }
    }
}