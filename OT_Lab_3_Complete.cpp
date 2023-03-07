#include<iostream>
#include<bits/stdc++.h>

using namespace std;

class Simplex{

    private:
        int rows, cols;
        //coefficients of all variables
        std::vector <std::vector<double> > A;
        //constants of constraints
        std::vector<double> B;
        //coefficients of objective function
        std::vector<double> C;

        std::vector<int> free;

        double maximum;

        int free_var;

        bool isUnbounded;

    public:
        Simplex(std::vector <std::vector<double> > matrix,std::vector<double> b ,std::vector<double> c, std::vector<int> f, int temp){
            maximum = 0;
            isUnbounded = false;
            rows = matrix.size();
            cols = matrix[0].size();
            A.resize( rows , vector<double>( cols , 0 ) );
            B.resize(b.size());
            C.resize(c.size());
            free.resize(f.size());
            free_var = temp;

            //pass A[][] values to matrix
            for(int i= 0;i<rows;i++){             
                for(int j= 0; j< cols;j++ ){
                    A[i][j] = matrix[i][j];
                }
            }

            for(int i = 0; i<f.size(); i++){
                free[i] = f[i];
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
                cout<<"-x-x-x-x-x- Error unbounded -x-x-x-x-x-\n"<<endl;
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
                print();
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

        //print current basis
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
            double minimum = INT16_MAX;
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
            int n= 0;

            cout<<"\n\nInitial array (Not optimal) :\n";
            print();
            cout<<"\n";
            //cout<<"Final array(Optimal solution) :\n";
            /*
            while(!end){
                bool result = simplexAlgorithmCalculataion();
                if(result==true){
                    end = true;
                    }
            }
            cout<<"Answers for Constraints of variables\n";
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
                    cout<<"x"<<index+1<<": "<<B[index]<<"\n"; 
                }
                else{
                    cout<<"x"<<index+1<<": "<<0<<"\n";
                }
            }
            cout<<"\nMaximum value: "<<maximum<<"\n";
            */
           while(!end){
                bool result = true;
                
                int pivot_column = findPivotColumn();
                int pivot_row = findPivotRow(pivot_column);
                if(pivot_column<cols && pivot_row<rows){
                    cout<<"Pivot is located at : ("<<pivot_row+1<<","<<pivot_column+1<<")\n\n";
                }
                cout<<"-----------------------------------------------------\n\n";

                if(!checkOptimality()){
                    result = simplexAlgorithmCalculataion();
                }
                n++;

                if(result==true){
                    end = true;
                    cout<<"Final array :\n";
                    print();
                    cout<<"Final values of variables :\n";
                }else{

                    cout<<"Array at iteration "<<n<<" :\n";
                    print();

                    cout<<"Values of variables after "<<n<<" iterations :\n";
                }

                int free_count = 0,temp = 0;
                for(int i=0;i< A[0].size(); i++){  
                    int count0 = 0;
                    int index = 0;
                    if(temp == 2){
                        temp = 0;
                        free_count++;
                    }
                    for(int j=0; j< rows; j++){
                        if(A[j][i]==0.0){
                            count0 += 1;
                        }
                        else if(A[j][i]==1){
                            index = j;
                        }
                    }
                    if(count0 == rows -1 && A[index][i] == 1){
                        if(free[i-free_count-temp] == 1 && i < cols-rows){
                            cout<<" BV x"<<i+1-free_count-temp<<temp+1<<": "<<B[index]<<"\n";
                            temp++;
                        }else{
                            cout<<" BV  x"<<i+1-free_count<<": "<<B[index]<<"\n";
                        }
                    }
                    else{
                        if(free[i-free_count-temp] == 1 && i < cols-rows){
                            cout<<"NBV x"<<i+1-free_count-temp<<temp+1<<": "<<0<<"\n";
                            temp++;
                        }else{
                            cout<<"NBV  x"<<i+1-free_count<<": "<<0<<"\n";
                        }
                    }
                }
                if(result == false)
                    cout<<"\nValue of Objective Function: "<<maximum<<"\n\n";
                else
                    cout<<"\nMaximum Value of Objective Function : "<<maximum<<"\n\n";
            }
        }
};

int main(){
    //no. of equations
    int eqns;
    //no. of variables
    int vars;
    cout<<"Enter no. of equations : ";
    cin>>eqns;
    cout<<"Enter no. of variables : ";
    cin>>vars;
    cout<<"Enter no. of free variables :";
    int free_var;
    cin>>free_var;
    vector<int> n(vars,0);
    if(free_var != 0){
        cout<<"\nEnter position of free variables (in single line) \n";
        cout<<"(eg. x1,x3>=0 & x2 free, input :  0 1 0)         :";
        for(int i = 0; i < vars; i++){
            cin>>n[i];
        }

    } 
    //coeffiecient vector
    vector<vector <double>> A(eqns,vector<double> (vars+eqns+free_var,0));
    //constraints vector
    vector<double> b(eqns,0);
    //objective function vector
    vector<double> c(eqns+vars+free_var,0);

    //input
    int free_count;
    cout<<"\nEnter A and b :"<<"\n\n";
    for(int i = 0; i < eqns; i++){
        free_count = 0;
        for(int j = 0; j < vars+free_var; j++){
            cout<<"Enter A[ "<<i+1-free_count<<" , "<<j+1-free_count<<" ] : ";
            if(n[j] == 0){
                cin>>A[i][j];
            }else{
                int temp;
                cin>>temp;
                A[i][j++] = temp;
                A[i][j] = -temp;
                free_count++;
            }
        }
        cout<<"Enter B[ "<<i+1<<" ] : ";
        cin>>b[i];
        cout<<"\n";
    }
    free_count = 0;
    cout<<"Enter C :\n";
    cout<<"(e.g. if z=x1+x2, enter C[0] : -1 C[1] : -1)\n\n";
    for(int i = 0; i < vars+free_var; i++){
        cout<<"Enter C[ "<<i+1-free_count<<" ] : ";
        if(n[i] == 1){
            int temp;
            cin>>temp;
            c[i++] = temp;
            c[i] = -temp;
            free_count++;
            continue;
        }
        cin>>c[i];
    }

    for(int i = 0; i+vars+free_var < vars+eqns+free_var; i++){
        A[i][i+vars+free_var] = 1;
    }
    
    cout<<"\nEntered function : \n";
    for(int i = 0; i < eqns + vars + free_var; i++){
        cout<<c[i]<<"\t";
    }

    cout<<"\n\nA : \n";
    for(int i = 0; i < eqns; i++){
        for(int j = 0; j < eqns+vars+free_var; j++){
            cout<<A[i][j]<<"\t";
        }
        cout<<"\n";
    }

    cout<<"\nB : \n";
    for(int i = 0; i < eqns; i++){
        cout<<b[i]<<"\t";
    }
    cout<<"\n";

    Simplex simplex(A,b,c,n,free_var);
    simplex.CalculateSimplex();
}