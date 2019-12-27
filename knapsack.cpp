#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <stdio.h>
#include <chrono>
#include <algorithm>
#include <limits>
#include <queue>
#include <thread>

struct item
{
    int id;         
    int profit;     
    int weight;     
    int bound;
    int level;
};


struct knapsack_problem
{
    std::vector<item> items;     
    int capacity;           
    int time_limit = 10000; 
};

struct bb_item {
    int maxProfit;     
    int maxWeight;     
    unsigned long duration;
};

unsigned long lastSolutionTime = 0;
//==============read_data
  knapsack_problem read_data(std::string& fn) {
    std::ifstream fin(fn); 
    if (!fin.good()) {     
        std::cerr << "Error opening file " << fn << std::endl;
        system("pause");
        exit(-1);         
    }

    knapsack_problem ks;  
                          
    int items_number;     
                          
    fin >> items_number;  
    for (int i = 0; i < items_number; i++) 
     
    {
        item an_item;     
        fin >> an_item.id;
        fin >> an_item.profit; 
        fin >> an_item.weight; 
        ks.items.push_back(an_item);
    }
    fin >> ks.capacity;         
    return ks; 
}

//=============print_knapsack_problem_info
  
void print_knapsack_problem_info(knapsack_problem& ks)
{
    std::cout << "Items=" << ks.items.size() << "   (id, Profit, Weight)\n";
    for (int i = 0; i < ks.items.size(); i++)            
                                                         
    {                                                    
        std::cout << ks.items[i].id << "  " << ks.items[i].profit << "  " << ks.items[i].weight << std::endl;
    }
    std::cout << "Capacity=" << ks.capacity << std::endl;
}

//================================================================
bool cmp(struct item a, struct item b)
{
    double r1 = (double)a.profit / a.weight;
    double r2 = (double)b.profit / b.weight;
    return r1 > r2;
}


 //============================================= greedy_approach
std::vector <item> greedy_approach(knapsack_problem& ks){
    auto start = std::chrono::high_resolution_clock::now();
    std::cout << "--- Greedy approach\n";
    sort(ks.items.begin(), ks.items.end(), cmp);

    int curWeight = 0; 
    double finalvalue = 0.0;   
    std::vector<item> solutionList;

    for (int i = 0; i < ks.items.size(); i++) {
        item tempItem = ks.items.at(i);
        if (curWeight + tempItem.weight <= ks.capacity)
        {
            solutionList.push_back(tempItem);
            curWeight += tempItem.weight;
            finalvalue += tempItem.profit;
        }

        else
        {
            int remain = ks.capacity - curWeight;
            finalvalue += tempItem.profit * ((double)remain / tempItem.weight);
            break;
        }
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
        lastSolutionTime = duration.count();
        if (lastSolutionTime > ks.time_limit) {
            std::cout << "Time out: t>" << ks.time_limit << "\n";
            return solutionList;
        }
    }

    std::cout << "maxProfit=" << finalvalue << "\n";
    std::cout << "maxWeight=" << curWeight << "\n";
    std::cout << "time=" << lastSolutionTime << "\n";
    return solutionList;
}


//================ Βrute Force Solver ======================
int get_profit(knapsack_problem& ks, std::vector<item>& sol){
    int total_profit = 0;
    int total_weight = 0;
    for (item an_item : sol)
    {
        total_weight += an_item.weight;
        if (total_weight > ks.capacity)
            return -1;
        total_profit += an_item.profit;
    }
    return total_profit;
}
//================================================================= brute_force_solver
std::vector <item> brute_force_solver(knapsack_problem& ks) {
    std::cout << "--- Brute force\n";
    auto start = std::chrono::high_resolution_clock::now();
    std::vector<item> result;
    int max_profit = -1;
    int n = ks.items.size();
    int total = 1 << n;
    if (n >= 31)  total = std::numeric_limits<int>::max();
    for (int i = 0;i < total;i++)
    {
        std::vector<item> sol;
        for (int j = 0;j < n;j++)
        {
            if ((i >> j) & 1)
                sol.push_back(ks.items[j]);

        }
        int profit = get_profit(ks, sol);
        if (profit > max_profit)
        {
            max_profit = profit;
            result = sol;

        }
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
        lastSolutionTime = duration.count();
        if (lastSolutionTime > 1000) {//ks.time_limit){
            std::cout << "Time out: t>" << ks.time_limit << "\n";
            return result;
        }
    }

    std::cout << "maxProfit=" << max_profit << "\n";
    int maxWeight = 0;
    for (auto myItem : result) maxWeight += myItem.weight;
    std::cout << "maxWeight=" << maxWeight << "\n";
    std::cout << "time=" << lastSolutionTime << "\n";

    return result;

}

//=========================
int bound(item u, int n, int W, item arr[]) {
    if (u.weight >= W) return 0;
    int profit_bound = u.profit;
    int j = u.level + 1;
    int totweight = u.weight;

    while ((j < n) && (totweight + arr[j].weight <= W)) {
        totweight += arr[j].weight;
        profit_bound += arr[j].profit;
        j++;
    }
    if (j < n) profit_bound += (W - totweight) * arr[j].profit / arr[j].weight;
    return profit_bound;
}

//===================branch_and_bound

bb_item branch_and_bound_action(int W, item arr[], int n, int timeLimit) {
    std::sort(arr, arr + n, cmp);
    std::queue <item> Q;
    item u, v;
    bb_item solutionItem;
    auto start = std::chrono::high_resolution_clock::now();
    u.level = -1;
    u.profit = u.weight = 0;
    Q.push(u);
    int maxProfit = 0;
    int maxWeight = 0;
    while (!Q.empty()) {
        u = Q.front();
        Q.pop();

        if (u.level == -1) v.level = 0;
        if (u.level == n - 1) continue;

        v.level = u.level + 1;

        v.weight = u.weight + arr[v.level].weight;
        v.profit = u.profit + arr[v.level].profit;
        v.id = arr[v.level].id;

        if (v.weight <= W && v.profit > maxProfit) {
            maxProfit = v.profit;
            maxWeight = v.weight;
        }

        v.bound = bound(v, n, W, arr);

        if (v.bound > maxProfit) {
            Q.push(v);
        }

        v.weight = u.weight;
        v.profit = u.profit;
        v.id = u.id;
        v.bound = bound(v, n, W, arr);
        if (v.bound > maxProfit) Q.push(v);

        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
        lastSolutionTime = duration.count();
        if (lastSolutionTime > timeLimit) {
            std::cout << "\nTime out: t>" << timeLimit << "\n";
            solutionItem.maxProfit = maxProfit;
            solutionItem.maxWeight = maxWeight;
            return solutionItem;
        }

    }
    solutionItem.maxProfit = maxProfit;
    solutionItem.maxWeight = maxWeight;
    std::cout << "maxProfit=" << maxProfit << "\n";
    std::cout << "maxWeight=" << maxWeight << "\n";
    std::cout << "time=" << lastSolutionTime << "\n";
    return solutionItem;
}
//=============================================================
bb_item branch_and_bound(knapsack_problem& ks) {

    bb_item solution;
    int n = ks.items.size();
    item* arr = new item[n];
    for (int i = 0;i < n;i++) arr[i] = ks.items.at(i);
    std::cout << "--- Branch and bound\n";
    solution = branch_and_bound_action(ks.capacity, arr, n, ks.time_limit);
   delete[] arr;

    return solution;
}

//========== Dynamic Algorithm ======================
int max(int x, int y) {return (x > y) ? x : y;}

std::vector <item> dynamicAlgorithm(knapsack_problem& ks) {
    std::cout << "--- Dynamic\n";
    auto start = std::chrono::high_resolution_clock::now();
    bool isTimeout = false;

    int W = ks.capacity;                      
    int n = ks.items.size();                  
    std::vector<int> w(n, 0);   
    std::vector<int> v(n, 0);   
    int m = 0;
    for (auto& myItem : ks.items) {
        w[m] = myItem.weight;
        v[m] = myItem.profit;
        m++;
    }

    int i, wt;
    int n_ = n + 1;
    int W_ = W + 1;

    std::vector<std::vector <int> > K(n_, std::vector<int>(W_, 0));
    for (i = 0; i <= n; i++) {                 

        for (wt = 0; wt <= W; wt++) {           
            if (i == 0 || wt == 0) K[i][wt] = 0; 
            else if (w[(i - 1)] <= wt)   
                K[i][wt] = max(v[(i - 1)] + K[i - 1][wt - w[(i - 1)]], K[i - 1][wt]);
            else
                K[i][wt] = K[i - 1][wt]; 
        }

        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
        lastSolutionTime = duration.count();
        if (lastSolutionTime > ks.time_limit) {
            std::cout << "Time out: t>" << ks.time_limit << "\n";
            isTimeout = true;
            break;
        }


    }
    lastSolutionTime = 0;
    std::cout << "end\n";
    // int bestProfit = K[n][W]        
    std::vector <item> solution;       
    if (!isTimeout) {
        while (n != 0) {               
                                         
            if (K[n][W] != K[n - 1][W]) { 
                                          
                item sol = ks.items.at((n - 1));
                solution.push_back(sol);        
                W = W - w[(n - 1)];             
                                                
            }
            n--;  
                  
        }

    }
    else {
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
        lastSolutionTime = duration.count();
    }




    for (int i = 0; i < n + 1; ++i) K[i].clear();
    K.clear();
    v.clear();
    w.clear();
    return solution;                          
}



//===============print_solution
void print_solution(std::vector<item> sol) {
    std::string output;
    int total_weight = 0;
    int total_profit = 0;
    char buffer[100];
    for (int i = 0;i < sol.size();i++) {
        std::cout << sol[i].id << " " << sol[i].profit << " " << sol[i].weight << "\n";
        total_weight += sol[i].weight;
        total_profit += sol[i].profit;
    }
    std::cout << "Total profit:" << total_profit << "\n";
    std::cout << "Total weight:" << total_weight << "\n";

}


//============export_solution
void export_solution(knapsack_problem& ks, std::vector<item> sol, std::string& fn){
    std::ofstream ofs(fn);
    if (!ofs.good()) {
        std::cerr << "Error accesing file" << fn << std::endl;
        exit(-1);
    }
    std::string output;
    int total_weight = 0;
    int total_profit = 0;
    for (int i = 0;i < sol.size();i++)
    {
        output += sol[i].id;
        output += " ";
        output += sol[i].profit;
        output += " ";
        output += sol[i].weight;
        output += " \n";

        ofs << output;
        total_weight += sol[i].weight;
        total_profit += sol[i].profit;
    }
    
    output = " ";
    output += total_profit;
    output += " ";
    output += total_weight;
    output += "\n";

    ofs << output;
    ofs.close();
}


//================= calculateProfitWeight
void calculateProfitWeight(std::vector<item>& solution, std::fstream& fout) {
    int w = 0;
    int v = 0;
    for (int i = 0; i < solution.size();i++) {
        item myItem = solution.at(i);
        w += myItem.weight;
        v += myItem.profit;
    }
    std::string s = ";" + std::to_string(w) + ";" + std::to_string(v);
    //char s[20];
    //sprintf(s, ";%d;%d\0", w, v);
    fout << s;
}


//=============create_CSV
 void create_CSV(int selection) {
    int n[] = { 10,50,100,500 };
    int r[] = { 50,100,500,1000 };
    int type[] = { 1,2,3,4 };
   // char filename[50];
    std::string filename_string;
    int fileCount = 0;

    std::string csv_path = "_results.csv";
    std::remove(csv_path.c_str());
    std::fstream fout;
    std::cout << "----- csv_path=" << csv_path << "\n";
    fout.open(csv_path, std::ios::out);
    if (!fout.good()) {
        std::cerr << "--- Error opening file " << csv_path << std::endl;   
        system("pause");
        exit(-1);     
    }


    std::string header = "INSTANCE";
    if ((selection == 0) || (selection == 1)) header.append(";WEIGHT_GR;VALUE_GR;EXECUTION_TIME_GR");
    if ((selection == 0) || (selection == 2)) header.append(";WEIGHT_BF;VALUE_BF;EXECUTION_TIME_BF");
    if ((selection == 0) || (selection == 3)) header.append(";WEIGHT_BB;VALUE_BB;EXECUTION_TIME_BB");
    if ((selection == 0) || (selection == 4)) header.append(";WEIGHT_DP;VALUE_DP;EXECUTION_TIME_DP");
    if ((selection == 0) || (selection == 5)) header.append(";WEIGHT_OT;VALUE_OT;EXECUTION_TIME_IP");
    if ((selection == 0) || (selection == 6)) header.append(";WEIGHT_IP;VALUE_IP;EXECUTION_TIME_OT");

    header.append("\n");
    fout << header;


    for (int k1 = 0;k1 < 4;k1++) {
        for (int k2 = 0;k2 < 4;k2++) {
            for (int k3 = 0;k3 < 4;k3++) {
                for (int k4 = 1;k4 <= 5;k4++) {
                    int n_ = n[k1];
                    int r_ = r[k2];
                    int type_ = type[k3];
                    int i = k4;
                    int S = 5;
                    std::vector<item> result;
                    
                    filename_string = "problem_" + std::to_string(n_) + "_" + std::to_string(r_) + "_" + std::to_string(type_) + "_" + std::to_string(i) + "_5.txt";
                   // std::string problems_full_path = "knapsack_problems/"+ filename_string;
                    fileCount++;
                    std::cout << "------------------" << fileCount << "-------------------\n";
                    std::cout << "filename= " << filename_string << "\n";
                    fout << ";" << filename_string;
                    knapsack_problem ks = read_data(filename_string);
                    std::cout << "### ks.capacity=" << ks.capacity << "\n";
                    std::cout << "### ks.size=" << ks.items.size() << "\n";

                  

                    //----  greedy_approach ---------
                    if ((selection == 0) || (selection == 1)) {
                        result.clear();
                        result = greedy_approach(ks);
                        calculateProfitWeight(result, fout);
                        fout << ";" << lastSolutionTime;
                    }

                    //----  Brute force ---------
                    if ((selection == 0) || (selection == 2)) {
                        result.clear();
                        result = brute_force_solver(ks);
                        calculateProfitWeight(result, fout);
                        fout << ";" << lastSolutionTime;
                    }

                    //----  Branch and bound ---------
                    if ((selection == 0) || (selection == 3)) {
                        result.clear();
                        bb_item solution = branch_and_bound(ks);
                        //calculateProfitWeight(result , fout);
                        fout << ";" << solution.maxWeight << ";" << solution.maxProfit;
                        fout << ";" << lastSolutionTime;

                    }

                    //----  Dynamic algorithm ---------
                    if ((selection == 0) || (selection == 4)) {
                        result.clear();
                        result = dynamicAlgorithm(ks);
                        calculateProfitWeight(result, fout);
                        fout << ";" << lastSolutionTime;
                    }

                    //----  Dynamic OR-Tools solver---------
                    if ((selection == 0) || (selection == 5)) {
                        fout.close();
                        std::string filenameWithArguments = "algorithm1.exe " + filename_string + " " + csv_path;
                        std::cout << "Filename: " << filenameWithArguments << "\n";
                        int retCode = system(filenameWithArguments.c_str());
                        fout.open(csv_path, std::ios::app);
                    }


                    //----  Integer OR-Tools solver ---------
                    if ((selection == 0) || (selection == 6)) {
                        fout.close();
                        std::string filenameWithArguments = "algorithm2.exe " + filename_string + " " + csv_path;
                        std::cout << "Filename: " << filenameWithArguments << "\n";
                        int retCode = system(filenameWithArguments.c_str());
                        fout.open(csv_path, std::ios::app);
                    }

                    fout << "\n";
                }
            }
        }
    }

    fout.close();
    std::cout << "Files: " << fileCount << "\n";
}

/*============================================================= main
  =============================================================
 */
int main(int argc, char* argv[]) {
    if (argc != 2) {
        std::cout << "Please enter one argument \n";
        exit(0);
    }

    int selection = atoi(argv[1]);
    
        
    if (selection == 0) {
        create_CSV(0);
        system("pause");
    }
    else if (selection == 1) {
        create_CSV(1);
        system("pause");
    }
    else if (selection == 2) {
        create_CSV(2);
        system("pause");
    }
    else if (selection == 3) {
        create_CSV(3);
        system("pause");
    }
    else if (selection == 4) {
        create_CSV(4);
        system("pause");
    }
    else if (selection == 5) {
        create_CSV(5);
        system("pause");
    }
    else if (selection == 6) {
        create_CSV(6);
        system("pause");
    }
    else {
        std::cout << "Enter an argument between 0 - 6\n";
        system("pause");
}




}