#include "FTSP.h"

void VNS_test_back()
{
    bool improve = true;
    vector<int> bin_sum(solution[0][0] + 1);
    for (int j = 1; j <= nj; j++)
    {
        bin_sum[solution[0][j]] += Job_Time[j];
    }
    int loop = 1;
    int count = 1;
    while (improve == true)
    {
        improve = false;
        if (improve == false)
        {
            for (int j = nj; j >= 1; j--)
            {
                for (int bin_2 = solution[0][0]; bin_2 >= 1; bin_2--)
                {
                    int bin_1 = solution[0][j];
                    int a = bin_sum[bin_1];
                    int b = bin_sum[bin_2];
                    if ((bin_1 != bin_2) && (due_date - Job_Time[j] >= b) && (b + Job_Time[j] > a))     
                    {
                        bool check_prece = true;                                                      
                        if (bin_1 < bin_2)
                        {
                            for (int l_id = 0; l_id < direct_sucs[j].size(); l_id++)
                            {
                                int l = direct_sucs[j][l_id];
                                if (solution[0][l] - bin_2 < 1)
                                {
                                    check_prece = false;
                                    break;
                                }
                            }
                        }
                        else
                        {
                            for (int l_id = 0; l_id < direct_preds[j].size(); l_id++)
                            {
                                int l = direct_preds[j][l_id];
                                if (bin_2 - solution[0][l] < 1)
                                {
                                    check_prece = false;
                                    break;
                                }
                            }
                        }
                        check_prece = check_prece && config_matr[solution[1][bin_2]][j];
                        if (check_prece)                                                       
                        {
                            solution[0][j] = bin_2;
                            bin_sum[bin_2] = bin_sum[bin_2] + Job_Time[j];
                            bin_sum[bin_1] = bin_sum[bin_1] - Job_Time[j];
                            improve = true;
                        }
                    }
                }
            }
        }
        if (improve == false)
        {
            for (int item_1 = nj; item_1 >= 1; item_1--)
            {
                for (int item_2 = item_1 - 1; item_2 >= 1; item_2--)
                {
                    int bin_1 = solution[0][item_1];
                    int bin_2 = solution[0][item_2];
                    if (bin_1 != bin_2)
                    {
                        int a = bin_sum[bin_1];
                        int b = bin_sum[bin_2];
                        int delta = Job_Time[item_1] - Job_Time[item_2];
                        if ((due_date >= a - delta) && (due_date >= b + delta))
                        {
                            bool check_prece = true;
                            if (bin_1 < bin_2)
                            {
                                for (int l_id = 0; (l_id < direct_sucs[item_1].size()) && (check_prece); l_id++)
                                {
                                    int l = direct_sucs[item_1][l_id];
                                    if (solution[0][l] - bin_2 < 1)
                                    {
                                        check_prece = false;
                                    }
                                }
                                for (int l_id = 0; (l_id < direct_preds[item_2].size()) && (check_prece); l_id++)
                                {
                                    int l = direct_preds[item_2][l_id];
                                    if (bin_1 - solution[0][l] < 1)
                                    {
                                        check_prece = false;
                                    }
                                }
                            }
                            else
                            {
                                for (int l_id = 0; (l_id < direct_sucs[item_2].size()) && (check_prece); l_id++)
                                {
                                    int l = direct_sucs[item_2][l_id];
                                    if (solution[0][l] - bin_1 < 1)
                                    {
                                        check_prece = false;
                                    }
                                }
                                for (int l_id = 0; (l_id < direct_preds[item_1].size()) && (check_prece); l_id++)
                                {
                                    int l = direct_preds[item_1][l_id];
                                    if (bin_2 - solution[0][l] < 1)
                                    {
                                        check_prece = false;
                                    }
                                }
                            }
                            check_prece = check_prece && config_matr[solution[1][bin_2]][item_1] && config_matr[solution[1][bin_1]][item_2];
                            if (check_prece == true)
                            {
                                solution[0][item_1] = bin_2;
                                solution[0][item_2] = bin_1;
                                bin_sum[bin_1] -= delta;
                                bin_sum[bin_2] += delta;
                                improve = true;
                            }
                        }
                    }
                }
            }

        }
        vector<int> empty_bin;
        empty_bin.resize(solution[0][0] + 1);
        empty_bin[0] = 1;
        for (int i = 1; i <= solution[0][0]; i++)
        {
            if (bin_sum[i] == 0)
            {
                empty_bin[empty_bin[0]] = i;
                empty_bin[0]++;
            }
        }
        empty_bin[0]--;
        count++;
        if (empty_bin[0] > 0)
        {
            solution[0][0] -= empty_bin[0];
            for (int j = 1; j <= nj; j++)
            {
                for (int id = 1; id <= empty_bin[0]; id++)
                {
                    if (id < empty_bin[0] && empty_bin[id] <= solution[0][j] && empty_bin[id + 1] > solution[0][j])
                    {
                        solution[0][j] -= id;
                    }
                    if (id == empty_bin[0] && empty_bin[id] <= solution[0][j])
                    {
                        solution[0][j] -= id;
                    }
                }
            }
            count = 0;
            bin_sum = vector<int>(solution[0][0] + 1);
            for (int j = 1; j <= nj; j++)
            {
                bin_sum[solution[0][j]] += Job_Time[j];
            }
            for (int emp = 1; emp <= empty_bin[0]; emp++) {

                solution[1].erase(solution[1].begin() + empty_bin[emp] - emp + 1);
            }
        }

        if (loop > 9999 || count > 1000) break;
        loop++;
        vector<int>().swap(empty_bin);
    }
}

void VNS_test_foward()
{
    bool improve = true;
    vector<int> bin_sum(solution[0][0] + 1);
    for (int j = 1; j <= nj; j++)
    {
        bin_sum[solution[0][j]] += Job_Time[j];
    }
    int loop = 1;
    int count = 1;
    while (improve == true)
    {
        improve = false;
        if (improve == false)
        {
            for (int j = 1; j <= nj; j++)
            {
                for (int bin_2 = 1; bin_2 <= solution[0][0]; bin_2++)
                {
                    int bin_1 = solution[0][j];
                    int a = bin_sum[bin_1];
                    int b = bin_sum[bin_2];
                    if ((bin_1 != bin_2) && (due_date - Job_Time[j] >= b) && (b + Job_Time[j] > a))   
                    {
                        bool check_prece = true;                                                      
                        if (bin_1 < bin_2)
                        {
                            for (int l_id = 0; l_id < direct_sucs[j].size(); l_id++)
                            {
                                int l = direct_sucs[j][l_id];
                                if (solution[0][l] - bin_2 < 1)
                                {
                                    check_prece = false;
                                    break;
                                }
                            }
                        }
                        else
                        {
                            for (int l_id = 0; l_id < direct_preds[j].size(); l_id++)
                            {
                                int l = direct_preds[j][l_id];
                                if (bin_2 - solution[0][l] < 1)
                                {
                                    check_prece = false;
                                    break;
                                }
                            }
                        }
                        check_prece = check_prece && config_matr[solution[1][bin_2]][j];
                        if (check_prece)                                                      
                        {
                            solution[0][j] = bin_2;
                            bin_sum[bin_2] = bin_sum[bin_2] + Job_Time[j];
                            bin_sum[bin_1] = bin_sum[bin_1] - Job_Time[j];
                            improve = true;
                        }
                    }
                }
            }
        }
        if (improve == false)
        {
            for (int item_1 = 1; item_1 <= nj; item_1++)
            {
                for (int item_2 = item_1 + 1; item_2 <= nj; item_2++)
                {
                    int bin_1 = solution[0][item_1];
                    int bin_2 = solution[0][item_2];
                    if (bin_1 != bin_2)
                    {
                        int a = bin_sum[bin_1];
                        int b = bin_sum[bin_2];
                        int delta = Job_Time[item_1] - Job_Time[item_2];
                        if ((due_date >= a - delta) && (due_date >= b + delta))
                        {
                            bool check_prece = true;
                            if (bin_1 < bin_2)
                            {
                                for (int l_id = 0; (l_id < direct_sucs[item_1].size()) && (check_prece); l_id++)
                                {
                                    int l = direct_sucs[item_1][l_id];
                                    if (solution[0][l] - bin_2 < 1)
                                    {
                                        check_prece = false;
                                    }
                                }
                                for (int l_id = 0; (l_id < direct_preds[item_2].size()) && (check_prece); l_id++)
                                {
                                    int l = direct_preds[item_2][l_id];
                                    if (bin_1 - solution[0][l] < 1)
                                    {
                                        check_prece = false;
                                    }
                                }
                            }
                            else
                            {
                                for (int l_id = 0; (l_id < direct_sucs[item_2].size()) && (check_prece); l_id++)
                                {
                                    int l = direct_sucs[item_2][l_id];
                                    if (solution[0][l] - bin_1 < 1)
                                    {
                                        check_prece = false;
                                    }
                                }
                                for (int l_id = 0; (l_id < direct_preds[item_1].size()) && (check_prece); l_id++)
                                {
                                    int l = direct_preds[item_1][l_id];
                                    if (bin_2 - solution[0][l] < 1)
                                    {
                                        check_prece = false;
                                    }
                                }
                            }
                            check_prece = check_prece && config_matr[solution[1][bin_2]][item_1] && config_matr[solution[1][bin_1]][item_2];
                            if (check_prece == true)
                            {
                                solution[0][item_1] = bin_2;
                                solution[0][item_2] = bin_1;
                                bin_sum[bin_1] -= delta;
                                bin_sum[bin_2] += delta;
                                improve = true;
                            }
                        }
                    }
                }
            }

        }
        vector<int> empty_bin;
        empty_bin.resize(solution[0][0] + 1);
        empty_bin[0] = 1;
        for (int i = 1; i <= solution[0][0]; i++)
        {
            if (bin_sum[i] == 0)
            {
                empty_bin[empty_bin[0]] = i;
                empty_bin[0]++;
            }
        }
        empty_bin[0]--;
        count++;
        if (empty_bin[0] > 0)
        {
            solution[0][0] -= empty_bin[0];
            for (int j = 1; j <= nj; j++)
            {
                for (int id = 1; id <= empty_bin[0]; id++)
                {
                    if (id < empty_bin[0] && empty_bin[id] <= solution[0][j] && empty_bin[id + 1] > solution[0][j])
                    {
                        solution[0][j] -= id;
                    }
                    if (id == empty_bin[0] && empty_bin[id] <= solution[0][j])
                    {
                        solution[0][j] -= id;
                    }
                }
            }
            count = 0;
            bin_sum = vector<int>(solution[0][0] + 1);
            for (int j = 1; j <= nj; j++)
            {
                bin_sum[solution[0][j]] += Job_Time[j];
            }
            for (int emp = 1; emp <= empty_bin[0]; emp++) {

                solution[1].erase(solution[1].begin() + empty_bin[emp] - emp + 1);
            }
        }

        if (loop > 9999 || count > 1000) break;
        loop++;
        vector<int>().swap(empty_bin);
    }
}

void VNS_flight()
{
    vector<vector<int>> test_scheduling(solution[0][0] + 1);
    for (int i = 1; i <= nj; i++)
    {
        test_scheduling[solution[0][i]].push_back(i);
    }
    vector<int> workload(num_config + 1);
    for (int j = 1; j <= solution[0][0]; j++)
    {
        workload[solution[1][j]]++;
    }
    bool if_move = false;
    if (if_move == false)
    {
        for (int j1 = 1; j1 <= solution[0][0]; j1++)
        {
            for (int j2 = j1 + 1; j2 <= solution[0][0]; j2++)
            {
                int config1 = solution[1][j1];
                int config2 = solution[1][j2];
                bool config_flag = true;
                for (int i = 0; (i < test_scheduling[j1].size()) && config_flag; i++)
                {
                    config_flag = config_flag && config_matr[config2][test_scheduling[j1][i]];
                }
                for (int i = 0; (i < test_scheduling[j2].size()) && config_flag; i++)
                {
                    config_flag = config_flag && config_matr[config1][test_scheduling[j2][i]];
                }
                if (config_flag)
                {
                    int tmp = solution[1][j1];
                    solution[1][j1] = solution[1][j2];
                    solution[1][j2] = tmp;
                }

            }
        }
    }
}

void Pert()
{
    srand(seed);
    bool change = false;
    int k_max = 8 > int(nj / 15) ? 8 : int(nj / 15);
    int rho_iter = k_max;

    vector<int> X_copy(solution[0]);
    vector<int> Y_copy(solution[1]);
    vector<int> data_workload(num_config + 1);
    for (int j = 1; j <= X_copy[0]; j++)
    {
        data_workload[Y_copy[j]]++;
    }
    int count = 0;
    while (change == false && count < 5)
    {
        change = true;
        vector<int> remove_list;
        for (int i = 0; i < rho_iter;)
        {
            bool flag = true;
            int r = 1 + rand() % nj;
            for (int j = 0; j < remove_list.size(); j++)
            {
                if (r == remove_list[j])
                {
                    flag = false;
                    break;
                }
            }
            if (flag == true)
            {
                remove_list.push_back(r);
                i++;
            }
        }
        for (int i = 0; i < remove_list.size(); i++)
        {
            X_copy[remove_list[i]] = 0;
        }
        vector<int> empty_bin(X_copy[0] + 1);
        empty_bin.resize(X_copy[0] + 1);
        for (int i = 1; i <= nj; i++)
        {
            empty_bin[X_copy[i]] = 1;
        }
        vector<int> empty_bin_2 = { 0 };
        for (int j = 1; j <= X_copy[0]; j++)
        {
            if (empty_bin[j] == 0)
            {
                data_workload[Y_copy[j]]--;
                empty_bin_2.push_back(j);
                empty_bin_2[0]++;
            }
        }
        X_copy[0] -= empty_bin_2[0];
        for (int j = 1; j <= nj; j++)
        {
            for (int id = 1; id <= empty_bin_2[0]; id++)
            {
                if (id < empty_bin_2[0] && empty_bin_2[id] <= X_copy[j] && empty_bin_2[id + 1] > X_copy[j])
                {
                    X_copy[j] -= id;
                }
                if (id == empty_bin_2[0] && empty_bin_2[id] <= X_copy[j])
                {
                    X_copy[j] -= id;
                }
            }
        }

        for (int emp = 1; emp <= empty_bin_2[0]; emp++)
        {
            Y_copy.erase(Y_copy.begin() + empty_bin_2[emp] - emp + 1);
        }
        vector<int> bin_sum(X_copy[0] + 1);
        for (int i = 1; i <= nj; i++)
        {
            if (X_copy[i] > 0)
            {
                bin_sum[X_copy[i]] += Job_Time[i];
            }
        }
        vector<int> X_copy_2(X_copy);
        vector<int> Y_copy_2(Y_copy);
        for (int i = 0; i < remove_list.size(); i++)
        {
            int item = remove_list[i];
            bool flag = false;
            for (int bins = 1; bins <= X_copy_2[0]; bins++)
            {
                bool flag_cap = false;
                bool flag_prec = true;
                bool flag_config = false;
                if (bin_sum[bins] + Job_Time[item] <= due_date)
                {
                    flag_cap = true;
                    for (int j = 0; j < direct_preds[item].size(); j++)
                    {
                        if (bins - X_copy_2[direct_preds[item][j]] < 1)
                        {
                            flag_prec = false;
                            break;
                        }
                    }
                    for (int j = 0; (j < direct_sucs[item].size()) && (flag_prec); j++)
                    {
                        if (X_copy_2[direct_sucs[item][j]] - bins < 1)
                        {
                            flag_prec = false;
                        }
                    }
                    if (flag_prec)
                    {
                        flag_config = config_matr[Y_copy_2[bins]][item];

                        if (!flag_config)
                        {
                            for (int k = 1; k <= num_config; k++)
                            {
                                if (config_matr[k][item] && data_workload[k] < C_k[k])
                                {
                                    bool flag_adj = true;
                                    for (int j = 1; (j <= nj) && flag_adj; j++)
                                    {
                                        if (X_copy_2[j] == bins) flag_adj = flag_adj && config_matr[k][j];
                                    }
                                    if (flag_adj == true)
                                    {
                                        flag_config = true;
                                        data_workload[Y_copy_2[bins]]--;
                                        Y_copy_2[bins] = k;
                                        data_workload[k]++;
                                        break;
                                    }
                                }
                            }
                        }
                    }
                }
                if (flag_cap && flag_prec && flag_config)
                {
                    X_copy_2[item] = bins;
                    bin_sum[bins] += Job_Time[item];
                    flag = true;
                    break;
                }
            }
            if (flag == false)
            {
                for (int location = 1; location <= X_copy_2[0] + 1; location++)
                {
                    flag = true;
                    for (int j = 0; j < direct_preds[item].size(); j++)
                    {
                        if (location <= X_copy_2[direct_preds[item][j]])
                        {
                            flag = false;
                            break;
                        }
                    }
                    for (int j = 0; (j < direct_sucs[item].size()) && (flag); j++)
                    {
                        if (location > X_copy_2[direct_sucs[item][j]])
                        {
                            flag = false;
                            break;
                        }
                    }
                    if (flag == true)
                    {
                        for (int k = 1; k <= num_config; k++)
                        {
                            flag = false;
                            if (config_matr[k][item] && data_workload[k] <= C_k[k])
                            {
                                flag = true;
                                Y_copy_2.insert(Y_copy_2.begin() + location, k);
                                data_workload[k]++;
                                bin_sum.insert(bin_sum.begin() + location, Job_Time[item]);
                                for (int i = 1; i <= nj; i++)
                                {
                                    if (X_copy_2[i] >= location) X_copy_2[i]++;
                                }
                                X_copy_2[item] = location;
                                X_copy_2[0]++;
                                break;
                            }
                        }
                        break;
                    }
                }
            }
            change = change && flag;
            if (flag)
            {
                X_copy = X_copy_2;
                Y_copy = Y_copy_2;
            }
            else
            {
                break;
            }
        }
        if (change == false)
        {
            X_copy = solution[0];
            Y_copy = solution[1];
            rho_iter = rho_iter - 1;
            if (rho_iter == 1)
            {
                rho_iter = k_max;
                count++;
            }
        }
        else
        {
            solution[0] = X_copy;
            solution[1] = Y_copy;
        }
    }
}

void ILS_solve()
{
    clock_t  start_time = clock();
    int count1 = 0;
    while (best_solution[0][0] > LB && optimalBDP == 0 && count1 < 5)
    {
        /*for section 6.5
        VNS_test_back();
        VNS_test_foward();
        VNS_flight();*/
        int count2 = 0;
        while (count2 < 10)
        {
            VNS_test_back();
            VNS_test_foward();
            VNS_flight();
            count2++;
        }
        if (best_solution[0][0] > solution[0][0]) best_solution = solution;
        Pert();
        count1++;
    }
    t_ILS = (double)(clock() - start_time) / CLOCKS_PER_SEC;
}
