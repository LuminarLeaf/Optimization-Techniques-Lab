#include<iostream>
#include<bits/stdc++.h>

#define iter 100

using namespace std;

void Gauss_Seidel(double A[], double b[], double x[], int m){
    // solution array and initial assumption
    double x_prev[m], x_temp[m];
    for(int i = 0; i < m; i++){
        x[i] = 0;
        x_prev[i] = 0;
        x_temp[i] = 0;
    }

    // iterative loop
    bool equal = false;
    for(int i = 0; i < iter; i++){
        
        // gauss seidel
        for(int p = 0; p < m; p++){
            x[p] = (b[p]/A[p*m + p]);
            for(int q = 0; q < m; q++){
                if(p == q)
                    continue;
                x[p] = x[p] - ((A[p*m + q]/A[p*m + p]) * x_temp[q]);
                x_temp[p] = x[p];
            }
        }

        // //print solution
        // cout<<"Iteration "<<(i+1)<<" :"<<endl;
        // for(int p = 0; p < m; p++){
        //     printf("x%d = %f    ", p + 1, x[p]);
        // }
        // cout<<endl;

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

int fact(int n){
    if(n == 0)
        return 1;
    int res = 1;
    for(int i = 2; i <= n; i++)
        res = res * i;
    return res;
}

int nCr(int n, int r){
    return fact(n) / (fact(r) * fact(n-r));
}

void solve(double A[], double b[], double *X, int m, int n, int mCn, int arr[], int curr_pos, int  vars_left, int &sol_no){
    if(vars_left == 0){
        int count = 0;
        double a[m*m], x[m];
        for(int i = 0; i < m; i++){
            for(int j = 0; j < m; j++){
                a[j*m + i] = A[arr[count] + j*n];
            }
            count ++;
        }
        Gauss_Seidel(a,b,x,m);
        count = 0;
        for(int i = 0; i < m; i++){
            X[sol_no*n + arr[count]] = x[i];
            count++;
        }
        sol_no++;
    }else{
        for(int i = curr_pos; i < n - vars_left + 1; i++){
            arr[m-vars_left] = i;
            solve(A,b,X,m,n,mCn,arr,i+1,vars_left-1,sol_no);
        }
    }
}

int main(){
    //no. of equations
    int m;
    //no. of variables
    int n;
    cout<<"Enter number of equations: ";
    cin>>m;
    cout<<"Enter number of variables: ";
    cin>>n;

    //matrix a
    cout<<"Enter A:\n";
    double A[m*n];
    for(int i = 0; i < m ; i++){
        for(int j = 0; j < n; j++){
            cout<<"A["<<(i + 1)<<", "<<(j + 1)<<"] =";
            cin>>A[i*n + j];
        }
    }

    // matrix b
    cout<<"\nEnter B:\n";
    double b[m];
    for(int i = 0; i < m; i++){
        cout<<"B["<<i+1<<"] :";
        cin>>b[i];
    }
    //mCn = no. of solutions
    int mCn = nCr(n,m);
    int arr[m], temp = 0;
    double X[mCn*n];

    for(int i = 0; i < mCn*n; i++){
        X[i] = 0;
    }

    solve(A,b,X,m,n,mCn,arr,0,m,temp);

    cout<<"\n\n";
    bool isPositive = true;
    temp = 0;
    for(int i = 0; i < mCn; i++){
        isPositive = true;
        for(int j = 0; j < n; j++){
            if(X[j + i*n] < 0 || isnan(X[j + i*n]))
                isPositive = false;
        }
        if(isPositive == true){
            temp ++;
        }
        cout<<"Solution "<<i+1<<" : "<<endl;
        for(int j = 0; j < n; j++){
            cout<<"   x"<<j+1<<" = "<<X[j + i*n]<<endl;
        }
    }
    if(temp == 0){
        cout<<"No Feasible Solution with all x's +ve"<<endl;
    }else{
        cout<<"No. of Basic Feasible Solutions : "<<temp<<endl;
    }
}