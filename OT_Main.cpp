#include<bits/stdc++.h>

bool BigM = false;  // true for Big M, false for simplex

using namespace std;

class Simplex{

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
        //Zj
        vector<double> Z;
        //Zj - Cj
        vector<double> ZC;
        //temporary storage
        vector<double> Ctemp;

        double maximum;

        //no of free variables
        int free_var_count;
        //position of free variables
        vector<int> free_var_pos;

        //type of equation (s, g, e)
        vector<char> eqns_type;

        //type of variable
        vector<char> var_type;

        //variable type string
        vector<string> var_type_str;

        //position of basic variable
        vector<int> bv_pos;

        bool isUnbounded;
        bool isOptimal = false;

        //problem type 0 - max | 1 - min
        int problem_type;

        Simplex(int m,int n,vector<vector<double> > matrix, vector<double> b, vector<double> c, vector<int> f, int temp, vector<char> tipe, int ptype){
            maximum = 0;
            isUnbounded = false;
            rows = matrix.size();
            cols = matrix[0].size();
            A.resize(rows , vector<double>(cols, 0));
            B.resize(b.size());
            C.resize(c.size());
            Ctemp.resize(c.size());
            free_var_pos.resize(f.size());
            free_var_count = temp;
            var_type.resize(tipe.size());
            var_type_str.resize(tipe.size());
            problem_type = ptype;
            bv_pos.resize(rows);
            Z.resize(c.size());
            ZC.resize(c.size());

            eqns = m;
            vars = n;

            for(int i = 0; i < rows; i++){
                for(int j = 0; j < cols; j++){
                    A[i][j] = matrix[i][j];
                }
            }

            for(int i = 0; i < B.size(); i++){
                B[i] = b[i];
            }

            for(int i = 0; i < C.size(); i++){
                C[i] = c[i];
                Ctemp[i] = c[i];
            }

            for(int i = 0; i < free_var_pos.size(); i++){
                free_var_pos[i] = f[i];
            }

            for(int i = 0; i < var_type.size(); i++){
                var_type[i] = tipe[i];
            }

            for(int i = 0; i < bv_pos.size(); i++){
                bv_pos[i] = 0;
            }
            
            int freecount = 0, scount = 0, acount = 0;
            for(int i = 0; i < var_type.size(); i++){
                if(var_type[i] == 'f'){
                    var_type_str[i] = "x"+to_string(i+1-freecount)+"a";
                    var_type_str[i+1] = "x"+to_string(i+1-freecount)+"b";
                    freecount++;
                    i++;
                }else if(var_type[i] == 'x'){
                    var_type_str[i] = "x"+to_string(i+1-freecount);
                }else if(var_type[i] == 's'){
                    var_type_str[i] = "s"+to_string(scount+1);
                    scount++;
                }else if(var_type[i] == 'a'){
                    var_type_str[i] = "a"+to_string(acount+1);
                    acount++;
                }
            }

            init_bv_pos();
            find_Zj();
            find_Zj_Cj();
        }

        bool simplexAlgorithmCalculation(){
            if(checkOptimality() == true){
                return true;
            }
            int pivotColumn = findPivotColumn();

            if(isUnbounded == true){
                cout<<"\n-x-x-x-x-x- Unbounded -x-x-x-x-x-\n"<<endl;
                return true;
            }

            int pivotRow = findPivotRow(pivotColumn);

            cout<<"\nPivot Column : "<<pivotColumn+1<<endl;
            cout<<"Pivot Row     : "<<pivotRow+1<<endl;
            cout<<"Pivot Element : "<<A[pivotRow][pivotColumn]<<endl;
            doPivotting(pivotRow, pivotColumn);

            return false;
        }

        //TODO: Integrate Problem Type
        bool checkOptimality(){
            isOptimal = false;
            int count = 0;
            for(int i = 0; i < ZC.size(); i++){
                double value = ZC[i];
                if(problem_type == 0){
                    if(value >= 0)
                        count++;
                }else if(problem_type == 1){
                    if(value <= 0)
                        count++;
                }
            }
            if(count == ZC.size()){
                isOptimal = true;
            }
            return isOptimal;
        }

        void doPivotting(int pivotRow, int pivotColumn){
            double pivotElement = A[pivotRow][pivotColumn];
            for(int i = 0; i < A[pivotRow].size(); i++){
                A[pivotRow][i] = A[pivotRow][i] / pivotElement;
            }
            B[pivotRow] = B[pivotRow] / pivotElement;

            for(int i = 0; i < A.size(); i++){
                if(i != pivotRow){
                    double temp = A[i][pivotColumn];
                    for(int j = 0; j < A[i].size(); j++){
                        A[i][j] = A[i][j] - (temp * A[pivotRow][j]);
                    }
                    B[i] = B[i] - (temp * B[pivotRow]);
                }
            }

            bv_pos[pivotRow] = pivotColumn;
            find_Zj();
            find_Zj_Cj();
        }

        void init_bv_pos(){
            int count0, count1, temp = 0;
            for(int i = vars+free_var_count; i < cols; i++){
                count1 = count0 = 0;
                for(int j = 0; j < rows; j++){
                    if(A[j][i] == 1){
                        count1 ++;
                    }else if(A[j][i] == 0){
                        count0 ++;
                    }
                }
                if(count1 == 1 && count0 == rows-1){
                    bv_pos[temp++] = i;
                }
            }
        }

        //TODO: Integrate Problem Type
        void find_Zj(){
            for(int i = 0; i < cols; i++){
                Z[i] = 0;
                for(int j = 0; j < rows; j++){
                    // if(problem_type == 0)
                        Z[i] += A[j][i] * C[bv_pos[j]];
                    // else if(problem_type == 1)
                    //     Z[i] += A[j][i] * C[bv_pos[j]] * -1;
                }
            }
        }

        //TODO: Integrate Problem Type
        void find_Zj_Cj(){
            for(int i = 0; i < cols; i++){
                // if(problem_type == 0)
                    ZC[i] = Z[i] - C[i];
                // else
                //     ZC[i] = C[i] - Z[i];
            }
        }

        void print(){
            cout<<"\nA : "<<endl;
            for(int i = 0; i < rows; i++){
                for(int j = 0; j < cols; j++){
                    cout<<A[i][j]<<"\t";
                }
                cout<<endl;
            }
            cout<<"\nB : "<<endl;
            for(int i = 0; i < B.size(); i++){
                cout<<B[i]<<"\t";
            }
            cout<<endl;
            cout<<"\nC : "<<endl;
            for(int i = 0; i < C.size(); i++){
                cout<<C[i]<<"\t";
            }
            cout<<endl;
        }

        void print_tableau(){
            int freecount = 0;
            cout<<"\tC : \t";
            for(int i = 0; i < C.size(); i++){
                cout<<C[i]<<"\t";
            }
            cout<<endl;
            cout<<"CBV\tBV\t";
            int scount = 0, acount = 0;
            for(int i = 0; i < var_type_str.size(); i++){
                cout<<var_type_str[i]<<"\t";
            }
            cout<<"Bx"<<endl;
            for(int i = 0; i < rows; i++){
                cout<<C[bv_pos[i]]<<"\t"<<var_type_str[bv_pos[i]]<<"\t";
                for(int j = 0; j < cols; j++){
                    cout<<precision(A[i][j], 4)<<"\t";
                }
                cout<<B[i]<<endl;
            }
            cout<<"\tZj :\t";
            for(int i = 0; i < Z.size(); i++){
                cout<<Z[i]<<"\t";
            }
            cout<<endl;
                cout<<"\tZj-Cj :\t";
            for(int i = 0; i < ZC.size(); i++){
                cout<<ZC[i]<<"\t";
            }
            cout<<endl;
        }

        double precision(double num, int precision){
            double temp = num * pow(10, precision);
            temp = round(temp);
            temp = temp / pow(10, precision);
            return temp;
        }

        void print_Solution(){
            int count1 = 0, count0 = 0, index = 0;
            double obj_value = 0;
            for(int i = 0; i < vars+free_var_count; i++){
                count1 = count0 = 0;
                for(int j = 0; j < rows; j++){
                    if(A[j][i] == 1){
                        count1 ++;
                        index = j;
                    }else if(A[j][i] == 0){
                        count0 ++;
                    }
                }
                if(ZC[i] == 0){
                    count0++;
                }
                if(count1 == 1 && count0 == rows){
                    cout<<"BV "<<var_type_str[i]<<"\t= "<<B[index]<<endl;
                    obj_value += C[bv_pos[index]] * B[index];
                }else{
                    cout<<"NBV "<<var_type_str[i]<<"\t= "<<0<<endl;
                }
            }
            cout<<"\nObjective Value : "<<obj_value<<endl;
        }

        //TODO: Integrate Problem Type
        int findPivotColumn(){
            int location = 0;
            double min0_max1 = ZC[0];

            for(int i = 1; i < ZC.size(); i++){
                if(problem_type == 0){
                    if(ZC[i] < min0_max1){
                        min0_max1 = ZC[i];
                        location = i;
                    }
                }else if(problem_type == 1){
                    if(ZC[i] > min0_max1){
                        min0_max1 = ZC[i];
                        location = i;
                    }
                }
            }
            return location;
        }

        int findPivotRow(int pivotColumn){
            double positiveValues[rows+1];
            vector<double> result(rows+1, 0);
            int negativeValueCount = 0;

            for(int i = 0; i < rows; i++){
                if(A[i][pivotColumn] > 0){
                    positiveValues[i] = A[i][pivotColumn];
                }else{
                    positiveValues[i] = 0;
                    negativeValueCount++;
                }
            }
            if(ZC[pivotColumn] < 0){
                negativeValueCount++;
            }

            if(negativeValueCount == rows){
                isUnbounded = true;
            }else{
                for(int i = 0; i < rows; i++){
                    double value = positiveValues[i];
                    if(value > 0){
                        result[i] = B[i] / value;
                    }else{
                        result[i] = 0;
                    }
                }
            }
            double minimum = INT_MAX;
            int location = 0;
            for(int i = 0; i < result.size(); i++){
                if(result[i] < minimum && result[i] > 0){
                    minimum = result[i];
                    location = i;
                }
            }
            return location;
        }

        void doIterations(){
            bool end = false;
            int iteration = 0;
            while(!end){
                bool result = true;
                if(!checkOptimality()){
                    result = simplexAlgorithmCalculation();
                }
                iteration++;
                cout<<"\n\n";
                if(result){
                    end = true;
                    cout<<"Final Tableau : "<<endl;
                    print_tableau();
                    cout<<"\nSolution : "<<endl;
                    print_Solution();
                }else{
                    cout<<"Iteration "<<iteration<<" : "<<endl;
                    print_tableau();
                    cout<<"\nValue of variables at "<<iteration<<" iterations : "<<endl;
                    print_Solution();
                }
                if(iteration > 20){
                    end = true;
                    cout<<"No solution found"<<endl;
                }
            }
        }

        void Solve(){
            bool end = false;
            cout<<"\n\nInitial Tableau : \n";
            print_tableau();
            cout<<"\n";
            doIterations();
        }
};


int main(){
    //no of equations
    int eqns;
    //no of variables
    int vars;
    //min/max
    int prob_type;
    cout<<"Enter 0 to Maximise, 1 to Minimise : ";
    cin>>prob_type;
    cout<<"Enter no. of equations : ";
    cin>>eqns;
    cout<<"Enter no. of variables : ";
    cin>>vars;

    int free_count;
    int temp = 0;

    cout<<"Enter no. of free variables : ";
    int free_var_count;
    cin>>free_var_count;
    vector<int> free_var_pos(vars,0);

    if(free_var_count != 0){
        cout<<"\nEnter position of free variables (in single line) \n";
        cout<<"(eg. x1,x3>=0 & x2 free, input :  0 1 0)         : ";
        for(int i = 0; i < vars; i++){
            cin>>free_var_pos[i];
        }
    }

    cout<<"\nEnter equation type (eg eq1 : <= eq2 : >= eq3 : =, input : s g e) : ";
    vector<char> eqns_type(eqns,0);
    int artificial_vars = 0;
    temp = 0;
    for(int i = 0; i < eqns; i++){
        cin>>eqns_type[i];
        if(eqns_type[i] != 's' && eqns_type[i] != 'g' && eqns_type[i] != 'e'){
            cout<<"Invalid input\n";
            return 0x00a;
        }
        if(eqns_type[i] == 'g'){
            artificial_vars++;
            temp++;
        }
        if(eqns_type[i] == 'e'){
            temp++;
        }
    }
    if(artificial_vars == eqns && prob_type == 0){
        cout<<"Unbounded as all equations are of the > type"<<endl;
        return 0x00a;
    }
    if(temp > 0){
        BigM = true;
        cout<<"Big M method used since number of artificial variables > 0"<<endl;
    }

    vector<char> var_type(vars+eqns+free_var_count+artificial_vars, '0');
    temp = 0;
    for(int i = 0; i < free_var_pos.size();i++){
        if(free_var_pos[i] == 1){
            var_type[i+temp] = 'f';
            var_type[i+temp +1] = 'f';
            temp++;
        }else{
            var_type[i+temp] = 'x';
        }
    }

    //coefficitent vector
    vector<vector<double> > A(eqns,vector<double> (vars+eqns+free_var_count+artificial_vars, 0));
    //constraints vector
    vector<double> b(eqns, 0);
    //objective function vector
    vector<double> c(eqns+vars+free_var_count+artificial_vars, 0);

    //input
    cout<<"Enter A and b : \n"<<endl;
    for(int i = 0; i < eqns; i++){
        free_count = 0;
        for(int j = 0; j < vars+free_var_count; j++){
            cout<<"Enter A[ "<<i+1<<" , "<<j+1-free_count<<" ] : ";
            if(free_var_pos[j - free_count] == 0){
                cin>>A[i][j];
            }else{
                int temp;
                cin>>temp;
                //break xi into a - b when free variable
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
    cout<<"(e.g. if z=x1+2(x2), enter C[0] : 1 C[1] : 2)\n\n";
    for(int i = 0; i < vars+free_var_count; i++){
        cout<<"Enter C[ "<<i+1-free_count<<" ] : ";
        int temp;
        cin>>temp;
        //TODO: Minimization
        // if(prob_type == 1){
        //     temp = -temp;
        // }
        if(free_var_pos[i] == 1){
            c[i++] = temp;
            c[i] = -temp;
            free_count++;
            continue;
        }
        c[i] = temp;
    }

    temp = 0;
    for(int i = 0; i+vars+free_var_count < vars+eqns+free_var_count+artificial_vars; i++){
        if(eqns_type[i-temp] == 's'){
            A[i-temp][i+vars+free_var_count] = 1;
            var_type[i+vars+free_var_count] = 's';
        }else if(eqns_type[i-temp] == 'g'){
            A[i-temp][i+vars+free_var_count] = -1;
            A[i-temp][i+vars+free_var_count+1] = 1;
            var_type[i+vars+free_var_count] = 's';
            var_type[i+vars+free_var_count+1] = 'a';
            if(BigM)
                c[i+vars+free_var_count+1] = 100000;
            i++;
            temp++;
        }else if(eqns_type[i-temp] == 'e'){
            A[i-temp][i+vars+free_var_count] = 1;
            var_type[i+vars+free_var_count] = 'a';
            if(BigM)
                c[i+vars+free_var_count] = 100000;
        }
    }

    Simplex simplex(eqns, vars, A, b, c, free_var_pos, free_var_count, var_type, prob_type);
    simplex.print();
    simplex.Solve();
}