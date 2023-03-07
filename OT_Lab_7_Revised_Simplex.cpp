#include<bits/stdc++.h>

bool BigM = false;  // true for Big M, false for simplex

using namespace std;

class RevisedSimplex{
    public:
        //no of equations
        int eqns;
        //no of variables
        int vars;
        //rows in A
        int rows;
        //columns in A
        int cols;
        //coefficients of variables
        vector<vector<double> > A;
        //Bx
        vector<double> B;
        //Cj
        vector<double> C;
        
        //Basis Vector
        vector<vector<double> > Basis;
        //Inverse of the basis vector
        vector<vector<double> > Basis_Inv;
        //Inverse of the new Basis Vector
        vector<vector<double> > Basis_New_Inv;
        //Cj values corresponding to basic variables
        vector<double> Cb;
        //Zj-Cj
        vector<double> ZjCj;
        //values of basic variables
        vector<double> Xb;
        //Y
        vector<double> Y;
        //alpha vector
        vector<double> Alpha;
        
        //value of objective function
        double Z;

        //position of basic variable
        vector<int> bv_pos;

        RevisedSimplex(vector<vector<double> > matrixA, vector<double> matrixB, vector<double> matrixC){
            rows = matrixA.size();
            cols = matrixA[0].size();

            eqns = rows;
            vars = cols-rows;
            
            A.resize(rows , vector<double>(cols, 0));
            B.resize(matrixB.size());
            C.resize(matrixC.size());
            bv_pos.resize(eqns);
            Basis.resize(eqns, vector<double> (eqns, 0) );
            Basis_Inv.resize(eqns, vector<double> (eqns, 0) );
            Basis_New_Inv.resize(eqns, vector<double> (eqns, 0) );
            Cb.resize(eqns);
            ZjCj.resize(eqns+vars);
            Xb.resize(eqns);
            Y.resize(eqns);
            Alpha.resize(eqns);

            for(int i = 0; i < rows; i++){
                for(int j = 0; j < cols; j++){
                    A[i][j] = matrixA[i][j];
                }
            }

            for(int i = 0; i < B.size(); i++){
                B[i] = matrixB[i];
            }

            for(int i = 0; i < C.size(); i++){
                C[i] = matrixC[i];
            }

            for(int i = 0; i < eqns; i++){
                Basis[i][i] = 1;
                Basis_Inv[i][i] = 1;
                Cb[i] = 0;
                Xb[i] = 0;
            }

            for(int i = 0; i < bv_pos.size(); i++){
                bv_pos[i] = 0;
            }

            for(int i = 0; i < eqns; i++){
                bv_pos[i] = vars+i;
            }
        }

        int search(vector<int> arr, int key){
            for(int i = 0; i < arr.size(); i++){
                if(arr[i] == key)
                    return true;
            }
            return false;
        }

        void calculateInverse(int LVI){
            //Ci
            vector<double> Ci(eqns, 0);
            for(int i = 0; i < eqns; i++){
                Ci[i] = Basis[i][LVI];
            }

            //e
            vector<double> E(eqns, 0);
            for(int i = 0; i < eqns; i++){
                for(int j = 0; j < eqns; j++){
                    E[i] += Basis_Inv[i][j] * Ci[j];
                }
            }

            //eta vector
            vector<double> eta(eqns, 0);
            for(int i = 0; i < eqns; i++){
                if(i == LVI){
                    eta[i] = 1 / E[i];
                }else{
                    eta[i] = -1 * E[i] / E[LVI];
                }
            }

            //Ei
            vector<vector<double> > Ei(eqns, vector<double> (eqns, 0));
            for(int i = 0; i < eqns; i++)
                Ei[i][i] = 1;
            for(int i = 0; i < eqns; i++){
                Ei[i][LVI] = eta[i];
            }

            //calculating new inverse
            for(int i = 0; i < eqns; i++){
                for(int j = 0; j < eqns; j++){
                    Basis_New_Inv[i][j] = 0;
                    for(int k = 0; k < eqns; k++){
                        Basis_New_Inv[i][j] += Ei[i][k] * Basis_Inv[k][j];
                    }
                }
            }
        }

        //print Basis, Basis_Inv, Cb, Xb, ZjCj
        void printMatrices(){
            cout<<"Basis Matrix : "<<endl;
            for(int i = 0; i < eqns; i++){
                for(int j = 0; j < eqns; j++){
                    cout<<Basis[i][j]<<"\t";
                }
                cout<<endl;
            }
            cout<<endl;

            cout<<"Inverse of Basis Matrix : "<<endl;
            for(int i = 0; i < eqns; i++){
                for(int j = 0; j < eqns; j++){
                    cout<<Basis_Inv[i][j]<<"\t";
                }
                cout<<endl;
            }
            cout<<endl;

            cout<<"Cb : "<<endl;
            for(int i = 0; i < eqns; i++){
                cout<<Cb[i]<<"\t";
            }
            cout<<endl<<endl;

            cout<<"Xb : "<<endl;
            for(int i = 0; i < eqns; i++){
                cout<<Xb[i]<<"\t";
            }
            cout<<endl<<endl;
        }

        void Solve(){
            cout<<"\n\nRevised Simplex Method :\n";
            cout<<"Initial Matrix : \n";
            printMatrices();
            bool end = false;
            int count_iter = 0;
            //Entering Variable Index
            int EVI;
            //Leaving Variable Index
            int LVI;
            do{
                //----------Step 1----------
                //calculating Y
                for(int i = 0; i < eqns; i++){
                    Y[i] = 0;
                    for(int j = 0; j < eqns; j++){
                        Y[i] += Cb[j] * Basis_Inv[j][i];
                    }
                }

                int minZC_index = -1;
                double minZC = INT_MAX;
                //calculating zj-cj
                for(int i = 0; i < eqns+vars; i++){
                    if(search(bv_pos, i)){
                        ZjCj[i] = 0;
                    }else{
                        int YP = 0;
                        for(int j = 0; j < eqns; j++){
                            YP += Y[j] * A[j][i];
                        }
                        ZjCj[i] = YP - C[i];
                    }
                    if(ZjCj[i] < minZC){
                        minZC = ZjCj[i];
                        minZC_index = i;
                    }
                }

                //checking if all Zj-Cj are positive
                end = true;
                for(int i = 0; i < eqns+vars; i++){
                    if(ZjCj[i] < 0)
                        end = false;
                }

                //calculating Xb
                for(int i = 0; i < eqns; i++){
                    Xb[i] = 0;
                    for(int j = 0; j < eqns; j++){
                        Xb[i] += Basis_Inv[i][j] * B[j];
                    }
                }

                //calculating Z (objective value)
                Z = 0;
                for(int i = 0; i < eqns; i++){
                    Z += Cb[i] * Xb[i];
                }

                //----------Step 2----------
                EVI = minZC_index;
                //only when not optimal
                if(end == false){
                    //calculating Alpha
                    for(int i = 0; i < eqns; i++){
                        Alpha[i] = 0;
                        for(int j = 0; j < eqns; j++){
                            Alpha[i] += Basis_Inv[i][j] * A[j][EVI];
                        }
                    }

                    int count = 0;
                    //checking for unboundedness
                    for(int i = 0; i < eqns; i++)
                        if(Alpha[i] < 0)
                            count++;
                    if(count == eqns){
                        cout<<"\n\n\n-x-x-x-x-x-Problem is Unbounded-x-x-x-x-x-";
                        exit(0);
                    }

                    int min_theta_index = -1;
                    double min_theta = INT_MAX;
                    //finding index of leaving variable
                    for(int i = 0; i < eqns; i++){
                        double theta = Xb[i] / Alpha[i];
                        if(theta < min_theta){
                            min_theta = theta;
                            min_theta_index = i;
                        }
                    }
                    //leaving variable is associated with the minimum theta
                    LVI = min_theta_index;

                    //updating the basis
                    for(int i = 0; i < eqns; i++){
                        Basis[i][LVI] = A[i][EVI];
                    }

                    //----------Step 3----------
                    // calculating inverse
                    calculateInverse(LVI);
                    
                    bv_pos[LVI] = EVI;

                    for(int i = 0; i < eqns; i++){
                        for(int j = 0; j < eqns; j++)
                            Basis_Inv[i][j] = Basis_New_Inv[i][j];
                        Cb[i] = C[bv_pos[i]];
                    }
                    count_iter++;
                    
                    cout<<"Value of Matrices after "<<count_iter<<" iterations : "<<endl;
                    printMatrices();
                    cout<<"Z = "<<Z<<endl;

                    if(count_iter == 50){
                        cout<<"\n\n\nToo many iterations. Exiting...\n\n";
                        exit(0);
                    }
                }else{
                    cout<<"\nFinal Value of Matrices : \n";
                    printMatrices();
                    cout<<"Maximum Value of Objective Function : "<<Z<<"\n\n";
                }
            }while(end == false);
        }
};


int main(){
    //no of equations
    int eqns;
    //no of variables
    int vars;
    cout<<"Enter no. of equations : ";
    cin>>eqns;
    cout<<"Enter no. of variables : ";
    cin>>vars;

    //coefficitent vector
    vector<vector<double> > A(eqns,vector<double> (vars+eqns, 0));
    //constraints vector
    vector<double> b(eqns, 0);
    //objective function vector
    vector<double> C(eqns+vars, 0);

    //input
    cout<<"Enter A and b : \n"<<endl;
    for(int i = 0; i < eqns; i++){
        for(int j = 0; j < vars; j++){
            cout<<"Enter A[ "<<i+1<<" , "<<j+1<<" ] : ";
                cin>>A[i][j];
        }
        cout<<"Enter B[ "<<i+1<<" ] : ";
        cin>>b[i];
        cout<<"\n";
    }

    for(int i = 0; i < eqns; i++){
        A[i][vars+i] = 1;
    }

    cout<<"\nEnter C :\n";
    cout<<"(e.g. if z=x1+2(x2), enter C[0] : 1 C[1] : 2)\n\n";
    for(int i = 0; i < vars; i++){
        cout<<"Enter C[ "<<i+1<<" ] : ";
        cin>>C[i];
    }

    cout<<"\nA : \n";
    for(int i = 0; i < eqns; i++){
        for(int j = 0; j < vars+eqns; j++){
            cout<<A[i][j]<<"\t";
        }
        cout<<"\n";
    }

    cout<<"\nb : \n";
    for(int i = 0; i < eqns; i++){
        cout<<b[i]<<"\t";
    }

    cout<<"\nC : \n";
    for(int i = 0; i < vars+eqns; i++){
        cout<<C[i]<<"\t";
    }
    cout<<"\n";


    RevisedSimplex simplex(A, b, C);
    simplex.Solve();
}