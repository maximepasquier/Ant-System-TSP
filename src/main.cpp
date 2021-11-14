#include <iostream>
#include <fstream>
#include <string>
#include <random>
#include <chrono>
#include <cmath>
#include <algorithm>
#include <iterator>

typedef struct Coords
{
    float X, Y;
} Coords;

typedef struct City_P
{
    int city_number;
    float city_probability;
} City_P;

int number_of_cities(std::string path);
void read_dat(std::string path, Coords cities_list[]);
void print_cities(Coords cities_list[], int cities_number);
float greedy_distance(Coords cities_list[], int cities_number);
float greedy_distance_shuffle(Coords cities_list[], int cities_number);
float dist(Coords current_city, Coords new_city_to_visit);
void AS(Coords cities_list[], int t_max, int m, int cities_number, float Q, float rho, float alpha, float beta);
float compute_p(Coords cities_list[], float **tao, int cities_number, int current_city, int iterator, bool already_visited[], float alpha, float beta);
int Roulette_method(std::vector<City_P> p);
template <typename M, typename S>
void print_matrix(int t, M matrix, S size);
float compute_tour_length(Coords cities_list[], int ant_tour[], int cities_number);
void print_solution(int solution_path[], float solution_path_distance, int cities_number);

int main()
{
    //* Get number of cities
    int cities_number = number_of_cities("./cities.dat");
    //* Init cities list
    Coords cities_list[cities_number];
    //* Get cities data
    read_dat("./cities.dat", cities_list);
    print_cities(cities_list, cities_number);

    //* Greedy minimal distance
    int greedy_iteration = 100;
    double Lnn = 0, Lnn_shuffle = 0;
    for (int i = 0; i < greedy_iteration; i++)
    {
        Lnn += greedy_distance(cities_list, cities_number);
        Lnn_shuffle += greedy_distance_shuffle(cities_list, cities_number);
    }
    Lnn /= greedy_iteration;
    Lnn_shuffle /= greedy_iteration;

    std::cout << "Lnn is : " << Lnn << std::endl;
    std::cout << "Lnn_shuffle is : " << Lnn_shuffle << std::endl;
    std::cout << std::endl;

    //* Ant System
    int t_max = 100;
    int m = 100;
    float Q = Lnn;
    float rho = 0.1;
    float alpha = 1;
    float beta = 5;

    //+ Tao matrix

    AS(cities_list, t_max, m, cities_number, Q, rho, alpha, beta);

    /*
    int ant_tour[cities_number] = {12, 8, 5, 7, 2, 0, 1, 13, 9, 3, 16, 11, 4, 17, 10, 14, 6, 15};
    float tour_length = compute_tour_length(cities_list, ant_tour, cities_number);
    std::cout << "tour_length is : " << tour_length << std::endl;
    */
}

void AS(Coords cities_list[], int t_max, int m, int cities_number, float Q, float rho, float alpha, float beta)
{
    //+ Tao matrix
    float **tao = new float *[cities_number];
    for (int i = 0; i < cities_number; i++)
    {
        tao[i] = new float[cities_number]();
    }

    //+ Set tao matrix to 1/Q
    for (size_t i = 0; i < cities_number; i++)
    {
        for (size_t j = 0; j < cities_number; j++)
        {
            tao[i][j] = 1 / Q;
        }
    }

    //* Solutions
    int solution_path[cities_number];
    float solution_path_distance = 1000;

    for (int t = 1; t <= t_max; t++)
    {
        for (int k = 1; k <= m; k++)
        {
            unsigned int seed = std::chrono::system_clock::now().time_since_epoch().count();
            std::default_random_engine generator(seed);
            std::uniform_int_distribution<int> distribution(0, cities_number - 1);
            //+ Define the starting point for the ant k
            int current_city = distribution(generator);
            //+ Declare the boolean array for visited cities of the ant k
            bool already_visited[cities_number];
            //+ Initialize the boolean array to false
            for (size_t i = 0; i < cities_number; i++)
            {
                already_visited[i] = false;
            }
            //+ Set origine city as visited
            already_visited[current_city] = true;
            //+ Declare ant tour
            int ant_tour[cities_number];
            //+ Initialize starting point
            ant_tour[0] = current_city;
            //+ Visit all the other cities
            for (size_t i = 1; i < cities_number; i++)
            {
                //+ Init vector of probabilities
                std::vector<City_P> p;
                for (int iterator = 0; iterator < cities_number; iterator++)
                {
                    if (!already_visited[iterator])
                    {
                        City_P probability = {iterator, compute_p(cities_list, tao, cities_number, current_city, iterator, already_visited, alpha, beta)};
                        p.push_back(probability);
                    }
                }
                //+ Choose a city destination according to the probabilities p
                int city_destination = Roulette_method(p);
                //+ Save tour of the ant
                ant_tour[i] = city_destination;
                //+ city_destination est donc visitée
                already_visited[city_destination] = true;
                //+ Update the current city
                current_city = city_destination;
            }
            //+ Compute length of the ant tour
            float tour_length = compute_tour_length(cities_list, ant_tour, cities_number);
            //+ Compute d_tao on the path of the ant
            for (int i = 0; i < cities_number - 1; i++)
            {
                int src_city = ant_tour[i];
                int dest_city = ant_tour[i + 1];
                tao[src_city][dest_city] += Q / tour_length;
            }
            //* Ant k has a better solution ?
            if (tour_length < solution_path_distance)
            {
                solution_path_distance = tour_length;
                std::copy(ant_tour, ant_tour + cities_number, solution_path);
            }
        }
        //+ Update tao matrix
        for (size_t i = 0; i < cities_number; i++)
        {
            for (size_t j = 0; j < cities_number; j++)
            {
                tao[i][j] = (1 - rho) * tao[i][j];
            }
        }
    }
    print_solution(solution_path, solution_path_distance, cities_number);
}

void print_solution(int solution_path[], float solution_path_distance, int cities_number)
{
    std::cout << "Best path is : ";
    for (size_t i = 0; i < cities_number; i++)
    {
        std::cout << solution_path[i] << " ";
    }
    std::cout << solution_path[0] << std::endl;
    std::cout << "Distance of this path is : " << solution_path_distance << std::endl;
}

float compute_tour_length(Coords cities_list[], int ant_tour[], int cities_number)
{
    float total_tour_length = 0;
    for (size_t i = 0; i < cities_number - 1; i++)
    {
        int src_city = ant_tour[i];
        int dest_city = ant_tour[i + 1];
        //std::cout << "dist entre : " << ant_tour[i] << " et " << ant_tour[i + 1] << " is : " << dist(cities_list[src_city], cities_list[dest_city]) << std::endl;
        total_tour_length += dist(cities_list[src_city], cities_list[dest_city]);
    }
    total_tour_length += dist(cities_list[ant_tour[0]], cities_list[ant_tour[cities_number - 1]]);
    //std::cout << "back to trace : " << dist(cities_list[ant_tour[0]], cities_list[ant_tour[cities_number - 1]]) << std::endl;
    return total_tour_length;
}

int Roulette_method(std::vector<City_P> p)
{
    unsigned int seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine generator(seed);
    std::uniform_real_distribution<double> distribution(0.0, 1.0);
    std::vector<City_P> probabilities_cumul;

    for (City_P item : p)
    {
        float tmp_proba;
        if (!probabilities_cumul.empty())
        {
            tmp_proba = probabilities_cumul.back().city_probability;
        }
        else
        {
            tmp_proba = 0;
        }
        probabilities_cumul.push_back({item.city_number, item.city_probability + tmp_proba});
    }
    double r = distribution(generator); // génération d'un nombre r
    for (City_P item : probabilities_cumul)
    {
        if (r < item.city_probability)
        {
            return item.city_number;
        }
    }
    std::cout << "ERROR  with roulette !" << std::endl;
    return -1;
}

float compute_p(Coords cities_list[], float **tao, int cities_number, int current_city, int iterator, bool already_visited[], float alpha, float beta)
{
    //+ Compute the sum for the normalization of p
    float normalizer = 0;
    //+ Parcourir les villes pour récupérer les non visitées
    for (size_t l = 0; l < cities_number; l++)
    {
        //+ Ville non visitée
        if (!already_visited[l])
        {
            //+ Sommer et additionner au normalisateur
            normalizer += pow(tao[current_city][l], alpha) * pow(1 / dist(cities_list[current_city], cities_list[l]), beta);
        }
    }
    return pow(tao[current_city][iterator], alpha) * pow(1 / dist(cities_list[current_city], cities_list[iterator]), beta) / normalizer;
}

void read_dat(std::string path, Coords cities_list[])
{
    std::ifstream cities_file;
    cities_file.open(path);
    std::string line;
    int iterator = 0;
    while (getline(cities_file, line))
    {
        std::string delimiter = " ";
        size_t position_delimiter = line.find(delimiter);
        line.erase(0, position_delimiter + delimiter.length());
        for (int i = 0; i < 2; i++)
        {
            size_t position_delimiter = line.find(delimiter);
            std::string key = line.substr(0, position_delimiter);
            switch (i)
            {
            case 0:
                cities_list[iterator].X = std::stof(key);
                break;
            case 1:
                cities_list[iterator].Y = std::stof(key);
                break;
            }
            line.erase(0, position_delimiter + delimiter.length());
        }
        iterator++;
    }
    cities_file.close();
}

template <typename M, typename S>
void print_matrix(int t, M matrix, S size)
{
    std::cout << std::endl;
    std::cout << "Matrix at iteration " << t << " is :" << std::endl;
    for (size_t i = 0; i < size; i++)
    {
        for (size_t j = 0; j < size; j++)
        {
            std::cout << matrix[i][j] << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

int number_of_cities(std::string path)
{
    std::ifstream cities;
    cities.open(path);
    std::string line;
    int cities_count = 0;
    while (getline(cities, line))
    {
        cities_count++;
    }
    cities.close();
    return cities_count;
}

void print_cities(Coords cities_list[], int cities_number)
{
    std::cout << "List of cities is : " << std::endl;
    for (int i = 0; i < cities_number; i++)
    {
        std::cout << "City number " << i << ", X = " << cities_list[i].X << ", Y = " << cities_list[i].Y << std::endl;
    }
    std::cout << std::endl;
}

float greedy_distance(Coords cities_list[], int cities_number)
{
    unsigned int seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine generator(seed);
    std::uniform_int_distribution<int> distribution(0, cities_number - 1);
    int r = distribution(generator);
    int origine = r;
    float total_distance = 0;
    bool already_visited[cities_number];
    for (int i = 0; i < cities_number; i++)
    {
        already_visited[i] = false;
    }
    already_visited[origine] = true;
    int current_city = origine;

    for (int i = 0; i < cities_number - 1; i++)
    {

        float min_distance = 1000;
        int min_city;
        for (int new_city_to_visit = 0; new_city_to_visit < cities_number; new_city_to_visit++)
            for (int a = 0; a < cities_number; a++)
            {
                if (!already_visited[new_city_to_visit])
                {
                    float local_distance = dist(cities_list[current_city], cities_list[new_city_to_visit]);
                    if (local_distance < min_distance)
                    {
                        min_distance = local_distance;
                        min_city = new_city_to_visit;
                    }
                }
            }
        //std::cout << "city number " << min_city << ", distance : " << min_distance << std::endl;
        total_distance += min_distance;
        current_city = min_city;
        already_visited[min_city] = true;
    }
    //+ Add le retourn sur le premier noeud
    total_distance += dist(cities_list[current_city], cities_list[origine]);
    return total_distance;
}

float greedy_distance_shuffle(Coords cities_list[], int cities_number)
{
    unsigned int seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine generator(seed);
    std::uniform_int_distribution<int> distribution(0, cities_number - 1);
    int r = distribution(generator);
    int origine = r;
    float total_distance = 0;
    bool already_visited[cities_number];
    for (int i = 0; i < cities_number; i++)
    {
        already_visited[i] = false;
    }
    already_visited[origine] = true;
    int current_city = origine;

    for (int i = 0; i < cities_number - 1; i++)
    {

        float min_distance = 1000;
        int min_city;
        int list_cities[cities_number];
        for (int b = 0; b < cities_number; b++)
        {
            list_cities[b] = b;
        }
        std::random_shuffle(&list_cities[0], &list_cities[cities_number]);
        /*
        for (size_t i = 0; i < cities_number; i++)
        {
            std::cout << list_cities[i] << std::endl;
        }
        std::cout << std::endl;
        */

        for (int iter = 0; iter < cities_number; iter++)
        {
            int new_city_to_visit = list_cities[iter];
            for (int a = 0; a < cities_number; a++)
            {
                if (!already_visited[new_city_to_visit])
                {
                    float local_distance = dist(cities_list[current_city], cities_list[new_city_to_visit]);
                    if (local_distance < min_distance)
                    {
                        min_distance = local_distance;
                        min_city = new_city_to_visit;
                    }
                }
            }
        }
        //std::cout << "city number " << min_city << ", distance : " << min_distance << std::endl;
        total_distance += min_distance;
        current_city = min_city;
        already_visited[min_city] = true;
    }
    //+ Add le retourn sur le premier noeud
    total_distance += dist(cities_list[current_city], cities_list[origine]);
    return total_distance;
}

float dist(Coords current_city, Coords new_city_to_visit)
{
    //std::cout << current_city.X << " " << new_city_to_visit.X << " : " << current_city.Y << " " << new_city_to_visit.Y << std::endl;
    float delta_X = current_city.X - new_city_to_visit.X;
    float delta_Y = current_city.Y - new_city_to_visit.Y;
    //std::cout << "delta_X : " << delta_X << ", delta_Y : " << delta_Y << std::endl;

    float distance = sqrt(pow(delta_X, 2) + pow(delta_Y, 2));

    //std::cout << "distance : " << distance << std::endl;
    //std::cout << std::endl;
    return distance;
}