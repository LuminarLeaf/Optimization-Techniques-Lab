#include<iostream>
#include<bits/stdc++.h>

using namespace std;

#define iter 100

class Simplex{

    private:
        int rows, cols;
        //coefficients of all variables
        std::vector <std::vector<double> > A;
        //constants of constraints
        std::vector<double> B;
        //coefficients of objective function
        std::vector<double> C;

        double maximum;

        bool isUnbounded;

    public:
        Simplex(std::vector <std::vector<double> > matrix,std::vector<double> b ,std::vector<double> c){
            maximum = 0;
            isUnbounded = false;
            rows = matrix.size();
            cols = matrix[0].size();
            A.resize( rows , vector<double>( cols , 0 ) );
            B.resize(b.size());
            C.resize(c.size());

            //pass A[][] values to matrix
            for(int i= 0;i<rows;i++){             
                for(int j= 0; j< cols;j++ ){
                    A[i][j] = matrix[i][j];
                }
            }

            //pass c[] values to B vector
            for(int i=0; i< c.size() ;i++ ){      
                C[i] = c[i] ;
            }

            //pass b[] values to B vector
            for(int i=0; i< b.size();i++ ){      
                B[i] = b[i];
            }
        }

        bool simplexAlgorithmCalculataion(){
            //check whether table is optimal,if optimal no need to process further
            if(checkOptimality()==true){
			    return true;
            }

            //find column which has pivot
            int pivotColumn = findPivotColumn();

            if(isUnbounded == true){
                cout<<"Error unbounded"<<endl;
			    return true;
            }

            //find row with pivot value
            int pivotRow = findPivotRow(pivotColumn);

            //form next table according to pivot value
            doPivotting(pivotRow,pivotColumn);

            return false;
        }

        bool checkOptimality(){
            bool isOptimal = false;
            int positivecount = 0;

            //check if coefficients of objective function are negative
            for(int i=0; i<C.size();i++){
                double value = C[i];
                if(value >= 0){
                    positivecount++;
                }
            }
            //if all constraints are positive table is optimal
            if(positivecount == C.size()){
                isOptimal = true;
            }
            return isOptimal;
        }

        void doPivotting(int pivotRow, int pivotColumn){
            //pivot value
            double pivotValue = A[pivotRow][pivotColumn];
            //column with pivot
            double pivotRowVals[cols];
            //row with pivot
            double pivotColVals[rows];
            //row after processing pivot value
            double rowNew[cols];

            //set maximum step by step
            maximum = maximum - (C[pivotColumn]*(B[pivotRow]/pivotValue));  
            //get row that has pivot value
            for(int i=0;i<cols;i++){
                pivotRowVals[i] = A[pivotRow][i];
            }
            //get column that has pivot value
            for(int j=0;j<rows;j++){
                pivotColVals[j] = A[j][pivotColumn];
            }

            //set row values that has pivot value divided by pivot value and put into new row
            for(int k=0;k<cols;k++){
                rowNew[k] = pivotRowVals[k]/pivotValue;
            }

            B[pivotRow] = B[pivotRow]/pivotValue;

            //process other coefficients in A array by subtracting
            for(int m=0;m<rows;m++){
                if(m !=pivotRow){
                    for(int p=0;p<cols;p++){
                        double multiplyValue = pivotColVals[m];
                        A[m][p] = A[m][p] - (multiplyValue*rowNew[p]);
                    }
                }
            }

            //process values of B array
            for(int i=0;i<B.size();i++){
                if(i != pivotRow){
                    double multiplyValue = pivotColVals[i];
                    B[i] = B[i] - (multiplyValue*B[pivotRow]);
                }
            }
            //least coefficient of constraints of objective function
            double multiplyValue = C[pivotColumn];
            for(int i=0;i<C.size();i++){
                C[i] = C[i] - (multiplyValue*rowNew[i]);
            }
            //replacing pivot row in new calculated A array
            for(int i=0;i<cols;i++){
                A[pivotRow][i] = rowNew[i];
            }
        }

        //print current A array
        void print(){
            for(int i=0; i<rows;i++){
                for(int j=0;j<cols;j++){
                    cout<<A[i][j] <<"\t";
                }
                cout<<B[i]<<"\n";
            }
            for(int i=0; i<cols; i++){
                cout<<C[i]<<"\t";
            }
            cout<<"\n\n";
        }

        //find least coefficients of constraints in objective function's position
        int findPivotColumn(){
            int location = 0;
            double minm = C[0];

            for(int i=1;i<C.size();i++){
                if(C[i]<minm){
                    minm = C[i];
                    location = i;
                }
            }
            return location;
        }

        //find row with pivot value
        int findPivotRow(int pivotColumn){
            double positiveValues[rows];
            std::vector<double> result(rows,0);
            int negativeValueCount = 0;

            for(int i=0;i<rows;i++){
                if(A[i][pivotColumn]>0){
                    positiveValues[i] = A[i][pivotColumn];
                }
                else{
                    positiveValues[i]=0;
                    negativeValueCount+=1;
                }
            }
            //checking unbound condition if all values are negative ones
            if(negativeValueCount==rows){
                isUnbounded = true;
            }
            else{
                for(int i=0;i<rows;i++){
                    double value = positiveValues[i];
                    if(value>0){
                        result[i] = B[i]/value;
                    }
                    else{
                        result[i] = 0;
                    }
                }
            }
            //find minimum's location of smallest item of B array
            double minimum = 999999999999;
            int location = 0;
            for(int i=0;i<sizeof(result)/sizeof(result[0]);i++){
                if(result[i]>0){
                    if(result[i]<minimum){
                        minimum = result[i];
                        location = i;
                    }
                }
            }
            return location;

        }

        void CalculateSimplex(){
            bool end = false;

            cout<<"\n\nInitial array (Not optimal) :\n";
            print();

            cout<<"\n";
            // cout<<"Final array(Optimal solution) :\n";

            // while(!end){
            //     bool result = simplexAlgorithmCalculataion();
            //     if(result==true){
            //         end = true;
            //         }
            // }
            // print();
            // cout<<"Final values of variables\n";

            cout<<"Simplex Table after 1 iteration :\n";
            simplexAlgorithmCalculataion();
            print();
            cout<<"Values of Variables after 1 iteraion :\n";

            //every basic column has values, get it form B array
            for(int i=0;i< A[0].size(); i++){  
                int count0 = 0;
                int index = 0;
                for(int j=0; j< rows; j++){
                    if(A[j][i]==0.0){
                        count0 += 1;
                    }
                    else if(A[j][i]==1){
                        index = j;
                    }
                }
                if(count0 == rows -1 ){
                    //every basic column has values, get it form B array
                    cout<<" BV x"<<i+1<<": "<<B[index]<<"\n"; 
                }
                else{
                    cout<<"NBV x"<<i+1<<": "<<0<<"\n";
                }
            }

            cout<<"\nValue of Objective Function: "<<maximum<<"\n";
        }

        void Variable_Values(){
            cout<<"Values of Variables:\n";

            //every basic column has values, get it form B array
            for(int i=0;i< A.size(); i++){  
                int count0 = 0;
                int index = 0;
                for(int j=0; j< rows; j++){
                    if(A[j][i]==0.0){
                        count0 += 1;
                    }
                    else if(A[j][i]==1){
                        index = j;
                    }
                }
                if(count0 == rows -1 ){
                    //every basic column has values, get it form B array
                    cout<<"x"<<index+1<<": "<<B[index]<<"\n"; 
                }
                else{
                    cout<<"x"<<index+1<<": "<<0<<"\n";
                }
            }
        }
};

void Gauss_Seidel(vector<vector<double>> A, vector<double> b, vector<double> &x, int m){
    vector<double> x_prev(m,0);
    vector<double> x_temp(m,0);
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
            x[p] = (b[p]/A[p][p]);
            for(int q = 0; q < m; q++){
                if(p == q)
                    continue;
                x[p] = x[p] - ((A[p][q]/A[p][p]) * x_temp[q]);
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

void solve(vector<vector<double>> A, vector<double> b, vector<vector<double>> &X, int m, int n, int mCn, vector<int> arr, int curr_pos, int  vars_left, int &sol_no){
    if(vars_left == 0){
        int count = 0;
        vector<vector<double>> a(m,vector<double> (m,0));
        vector<double> x(m,0);
        for(int i = 0; i < m; i++){
            for(int j = 0; j < m; j++){
                a[i][j] = A[j][arr[count]];
                //a[j][i] = A[arr[count]][j];
            }
            count ++;
        }
        Gauss_Seidel(a,b,x,m);
        count = 0;
        for(int i = 0; i < m; i++){
            X[sol_no][arr[count]] = x[i];
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
    int eqns;
    //no. of variables
    int vars;
    cout<<"Enter no. of equations : ";
    cin>>eqns;
    cout<<"Enter no. of variables : ";
    cin>>vars;

    //coeffiecient vector
    vector<vector <double> > A(eqns,vector<double> (vars+eqns,0));
    //constraints vector
    vector<double> b(eqns,0);
    //objective function vector
    vector<double> c(eqns+vars,0);

    //input
    cout<<"Enter A and b :"<<"\n\n";
    for(int i = 0; i < eqns; i++){
        for(int j = 0; j < vars; j++){
            cout<<"Enter A[ "<<i+1<<" , "<<j+1<<" ] : ";
            cin>>A[i][j];
        }
        cout<<"Enter B[ "<<i+1<<" ] : ";
        cin>>b[i];
        cout<<"\n";
    }
    cout<<"Enter C :\n";
    cout<<"(e.g. if z=x1+x2, enter C[0] : -1 C[1] : -1)\n\n";
    for(int i = 0; i < vars; i++){
        cout<<"Enter C[ "<<i+1<<" ] : ";
        cin>>c[i];
    }

    for(int i = 0; i+vars < vars+eqns; i++){
        A[i][i+vars] = 1;
    }
    
    cout<<"\nEntered function : \n";
    for(int i = 0; i < eqns + vars; i++){
        cout<<c[i]<<"\t";
    }

    cout<<"\n\nA : \n";
    for(int i = 0; i < eqns; i++){
        for(int j = 0; j < eqns+vars; j++){
            cout<<A[i][j]<<"\t";
        }
        cout<<"\n";
    }

    cout<<"\nB : \n";
    for(int i = 0; i < eqns; i++){
        cout<<b[i]<<"\t";
    }
    cout<<"\n\n";
    cout<<"Choices:\n1) List all BFS\n2) Initial Simplex table\n3) List of all variables in initial table\n4) Simplex table of the 1st Iteration\n\n";
    tell:
    cout<<"Enter choice : ";
    int choice;
    cin>>choice;
    Simplex simplex(A,b,c);
    switch (choice)
    {
    case 1:
        {
            int mCn = nCr(vars+eqns, eqns);
            int temp = 0;
            vector<int> arr(eqns,0);
            vector<vector<double>> X(mCn,vector<double> (vars+eqns,0));
            solve(A,b,X,eqns,vars+eqns,mCn,arr,0,eqns,temp);

            cout<<"\n\n";
            bool isPositive = true;
            temp = 0;
            for(int i = 0; i < mCn; i++){
                isPositive = true;
                for(int j = 0; j < vars+eqns; j++){
                    if(X[i][j] < 0 || isnan(X[i][j]))
                        isPositive = false;
                }
                if(isPositive == true){
                    temp ++;
                }
                cout<<"Solution "<<i+1<<" : "<<endl;
                for(int j = 0; j < vars+eqns; j++){
                    cout<<"   x"<<j+1<<" = "<<X[i][j]<<endl;
                }
            }
            if(temp == 0){
                cout<<"No Feasible Solution with all x's +ve"<<endl;
            }else{
                cout<<"No. of Basic Feasible Solutions : "<<temp<<endl;
            }
        }
        break;
    case 2:
        simplex.print();
        break;
    case 3:
        simplex.Variable_Values();
        break;
    case 4:
        simplex.CalculateSimplex();
        break;
    default:
        cout<<"Please enter correct choice.\n\n";
        goto tell;
        break;
    }
}