#include <algorithm>
#include <array>
#include <fstream>
#include <iostream>
#include <random>
#include <vector>

// INFO: Constants
#define GRID_SIZE 9
#define POPULATION_SIZE 1000
#define BEST_SIZE 200
#define RANDOM_SIZE 200
#define MUTATION_RATE 0.1
#define GENERATIONS_SIZE 10000
#define STAGNATION_THRESHOLD 100

// INFO: Types definitions
// Single element<value, predefined>
typedef std::pair<int, int16_t> elem;
// Single row
typedef std::vector<elem> gen;
// Table
typedef std::vector<gen> chrom;
// Vector of tables
typedef std::vector<chrom> population;

// INFO: Declarations
/*
 * @brief Read file with initial sudoku
 * @return initial sudoku in 2-D vector format
 */
chrom read_file();

/*
 * @brief Generate undefined parts of the chromosome (single row)
 * @param initial gene (separate row)
 * @return initial gene
 */
gen generate_gen(const gen& initial);

/*
 * @brief Generate undefined parts of the chromosome
 * @param initial chromosome (separate table)
 * @return Initial chromosome
 */
chrom generate_chromosome(const chrom& initial);

/*
 * @brief Generate initial population (set of tables)
 * @param initial sudoku
 * @return vector of initial chromosomes
 * */
population generate_init_population(int size, const chrom& initial);

/*
 * @brief Calculate fitness function value of generated solution
 * @param current solution (chromosome/table)
 * @return fitness function value
 */
double fitness(const chrom& current);

/*
 * @brief Uniform-Order crossover (swap elements of two chromosomes)
 * @param initial chromosome
 * @return crossovered chromosome
 */
std::array<chrom, 2> crossover(chrom& first, chrom& second);

/*
 * @brief Mutate given chromosome (table)
 * @param chromosome (table)
 */
void mutate(chrom& current);

/*
 * @brief Generate new population based of fitness function
 * @param current population
 * @return new population (next step of solutions)
 */
population selection(population& current);

/*
 * @brief Generate new population based crossover and mutation
 * @param current population
 * @return new population (next step of solutions)
 */
population generate_population(population current);

/*
 * @brief Print result in default output
 * @param current solution (with maximum fitness function value)
 */
void print_result(const chrom& current);

/*
 * Generate solution
 * @brief general function to collect all functionality
 * @return chromosome (either solution or the best one)
 */
chrom generate_solution();

chrom find_best(population& pop) {
  chrom best = pop[0];
  double best_fitness = fitness(best);
  for (const chrom& chrom : pop) {
    double fit = fitness(chrom);
    if (fit < best_fitness) {
      best = chrom;
      best_fitness = fit;
    }
  }
  return best;
}

population generate_population_from_solution(int size,
                                             const chrom& best_solution) {
  population pop;
  for (int i = 0; i < size; ++i) {
    chrom new_chrom = best_solution;
    mutate(new_chrom);
    pop.push_back(new_chrom);
  }
  return pop;
}

// INFO: Driver function
int main() {
  std::ios::sync_with_stdio(false);
  std::cin.tie(nullptr);

  chrom solution = generate_solution();
  print_result(solution);

  return 0;
}

// INFO: Definitions
chrom read_file() {
  std::fstream in("input.txt");
  chrom table(GRID_SIZE, gen(GRID_SIZE));

  char symb;
  for (int row = 0; row != GRID_SIZE; ++row) {
    for (int col = 0; col != GRID_SIZE; ++col) {
      in >> symb;
      if (symb != '-') {
        table[row][col] = {symb - '0', 1};
      } else {
        table[row][col] = {0, 0};
      }
    }
  }
  in.close();

  return table;
}

gen generate_gen(const gen& initial) {
  // Generate random row of elements
  gen new_gen = {{1, 0}, {2, 0}, {3, 0}, {4, 0}, {5, 0},
                 {6, 0}, {7, 0}, {8, 0}, {9, 0}};
  // Shuffle them to change order
  std::shuffle(new_gen.begin(), new_gen.end(),
               std::mt19937(std::random_device()()));

  // Set predefined elements back
  for (int i = 0; i != GRID_SIZE; ++i) {
    if (initial[i].second) {
      new_gen[i] = initial[i];
    }
  }
  return new_gen;
}

chrom generate_chromosome(const chrom& initial) {
  chrom new_chrom;
  for (int row = 0; row != GRID_SIZE; ++row) {
    new_chrom.push_back(generate_gen(initial[row]));
  }
  return new_chrom;
}

population generate_init_population(int size, const chrom& initial) {
  population population;
  for (int i = 0; i != size; ++i) {
    population.push_back(generate_chromosome(initial));
  }
  return population;
}

double fitness(const chrom& chrom) {
  double total_fitness = 0.0;

  for (int row = 0; row < GRID_SIZE; ++row) {
    int predefined_count = 0;
    for (int col = 0; col < GRID_SIZE; ++col) {
      if (chrom[row][col].second) {
        ++predefined_count;
      }
    }

    double weight = (9 - predefined_count) / 9.0;
    uint32_t bitmask = 0;
    for (int col = 0; col < GRID_SIZE; ++col) {
      int num = chrom[row][col].first - 1;
      if (bitmask & (1 << num)) {
        total_fitness += weight;
      } else {
        bitmask |= (1 << num);
      }
    }
  }

  for (int col = 0; col < 9; ++col) {
    int predefined_count = 0;
    for (int row = 0; row < 9; ++row) {
      if (chrom[row][col].second) {
        ++predefined_count;
      }
    }
    double weight = (9 - predefined_count) / 9.0;
    uint32_t bitmask = 0;
    for (int row = 0; row < 9; ++row) {
      int num = chrom[row][col].first - 1;  // 0 to 8
      if (bitmask & (1 << num)) {
        total_fitness += weight;
      } else {
        bitmask |= (1 << num);
      }
    }
  }

  for (int block_row = 0; block_row < 3; ++block_row) {
    for (int block_col = 0; block_col < 3; ++block_col) {
      int predefined_count = 0;
      for (int row = block_row * 3; row < (block_row + 1) * 3; ++row) {
        for (int col = block_col * 3; col < (block_col + 1) * 3; ++col) {
          if (chrom[row][col].second) {
            ++predefined_count;
          }
        }
      }
      double weight = (9 - predefined_count) / 9.0;
      uint32_t bitmask = 0;
      for (int row = block_row * 3; row < (block_row + 1) * 3; ++row) {
        for (int col = block_col * 3; col < (block_col + 1) * 3; ++col) {
          int num = chrom[row][col].first - 1;
          if (bitmask & (1 << num)) {
            total_fitness += weight;
          } else {
            bitmask |= (1 << num);
          }
        }
      }
    }
  }

  return total_fitness;
}

std::array<chrom, 2> crossover(chrom& first_parent, chrom& second_parent) {
  std::random_device dev;
  std::mt19937 rng(dev());
  std::uniform_int_distribution<> select(0, GRID_SIZE - 1);

  chrom first_child = first_parent;
  chrom second_child = second_parent;

  for (int row = 0; row != GRID_SIZE; ++row) {
    int positions[2]{select(rng), select(rng)};
    if (positions[1] < positions[0]) {
      std::swap(positions[0], positions[1]);
    }

    for (int col = 0; col != positions[0]; ++col) {
      if (!first_child[row][col].second && second_parent[row][col].second) {
        std::swap(first_child[row][col], second_parent[row][col]);
      }
    }
    for (int col = positions[0]; col != positions[1]; ++col) {
      if (!second_child[row][col].second && !first_parent[row][col].second) {
        std::swap(first_parent[row][col], second_child[row][col]);
      }
    }
    for (int col = positions[1]; col != GRID_SIZE; ++col) {
      if (!second_parent[row][col].second && !first_child[row][col].second) {
        std::swap(first_child[row][col], second_parent[row][col]);
      }
    }
  }

  std::array<chrom, 2> children = {first_child, second_child};
  return children;
}

void mutate(chrom& current) {
  std::random_device dev;
  std::mt19937 rng(dev());
  std::uniform_int_distribution<> select_block(0, 8);

  int swap_coords[4]{select_block(rng), select_block(rng), select_block(rng),
                     select_block(rng)};

  if (!current[swap_coords[0]][swap_coords[1]].second &&
      !current[swap_coords[2]][swap_coords[3]].second) {
    std::swap(current[swap_coords[0]][swap_coords[1]],
              current[swap_coords[2]][swap_coords[3]]);
  }
}

population tournament_selection(population& current) {
  std::sort(current.begin(), current.end(), [](const chrom& a, const chrom& b) {
    return fitness(a) < fitness(b);
  });

  population new_pool;

  std::random_device dev;
  std::mt19937 rng(dev());
  std::uniform_int_distribution<> select_general(0, POPULATION_SIZE);
  std::uniform_real_distribution<> select_prob(0.0, 1.0);

  for (int i = 0; i != BEST_SIZE; ++i) {
    new_pool.push_back(current[i]);
  }

  while (new_pool.size() < BEST_SIZE + RANDOM_SIZE) {
    int first_parent = select_general(rng) % new_pool.size();
    int second_parent = select_general(rng) % new_pool.size();

    auto children = crossover(new_pool[first_parent], new_pool[second_parent]);

    if (select_prob(rng) < MUTATION_RATE) {
      mutate(children[0]);
    }
    if (select_prob(rng) < MUTATION_RATE) {
      mutate(children[1]);
    }

    new_pool.push_back(children[0]);
    new_pool.push_back(children[1]);

    new_pool.push_back(children[0]);

    if (new_pool.size() < BEST_SIZE + RANDOM_SIZE) {
      new_pool.push_back(children[1]);
    }
  }

  return new_pool;
}

chrom generate_solution() {
  chrom initial = read_file();
  population chrom_population =
      generate_init_population(POPULATION_SIZE, initial);

  chrom global_best = find_best(chrom_population);
  double global_best_fitness = fitness(global_best);

  int stagnation_counter = 0;

  for (int i = 0; i < GENERATIONS_SIZE; ++i) {
    population next_population = tournament_selection(chrom_population);

    chrom new_best = find_best(next_population);
    double new_best_fitness = fitness(new_best);

    if (new_best_fitness < global_best_fitness) {
      global_best = new_best;
      global_best_fitness = new_best_fitness;
      stagnation_counter = 0;
    } else {
      stagnation_counter++;
      if (stagnation_counter >= STAGNATION_THRESHOLD) {
        // Restart population with global_best as base
        chrom_population =
            generate_population_from_solution(POPULATION_SIZE, global_best);
        stagnation_counter = 0;
      }
    }

    if (global_best_fitness < fitness(next_population.back())) {
      next_population.back() = global_best;
    }

    chrom_population = next_population;

    std::cout << "Generation: " << i
              << "\tBest fitness: " << global_best_fitness << "\n";
    if (global_best_fitness == 0) {
      std::cout << "Solution found!\n";
      break;
    }
  }

  return global_best;
}
void print_result(const chrom& current) {
  std::ofstream out("output.txt");

  if (out.is_open()) {
    for (int row = 0; row != GRID_SIZE; ++row) {
      for (int col = 0; col != GRID_SIZE; ++col) {
        int current_el = +current[row][col].first;
        out << current_el << " ";
      }
      out << "\n";
    }
  }

  out.close();
}
