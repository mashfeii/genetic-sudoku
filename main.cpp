#include <algorithm>
#include <fstream>
#include <iostream>
#include <numeric>
#include <random>
#include <vector>

#define GRID_SIZE 9
#define POPULATION_SIZE 100
#define MUTATION_RATE 0.1
#define CROSSOVER_RATE 0.9
#define MAX_POPULATIONS 100
#define LINE_PRODUCT 362880

// Single element<value, predefined>
typedef std::pair<unsigned char, unsigned char> elem;
// Single row
typedef std::vector<elem> gen;
// Table
typedef std::vector<gen> chrom;
// Vector of tables
typedef std::vector<chrom> population;

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
gen generate_gen(gen& initial);

/*
 * @brief Generate undefined parts of the chromosome
 * @param initial chromosome (separate table)
 * @return Initial chromosome
 */
chrom generate_chromosome(chrom& initial);

/*
 * @brief Generate initial population (set of tables)
 * @param initial sudoku
 * @return vector of initial chromosomes
 * */
population generate_init_population(int size, chrom& initial);

/*
 * @brief Calculate fitness function value of generated solution
 * @param current solution (chromosome/table)
 * @return fitness function value
 */
int fitness(chrom& current);

/*
 * @brief Uniform-Order crossover (swap elements of two chromosomes)
 * @param initial chromosome
 * @return crossovered chromosome
 */
chrom crossover(chrom& first, chrom& second);

/*
 * @brief Mutation
 * @param initial chromosome
 * @return mutated chromosome
 */
chrom mutate(chrom& initial);

int main() {
  std::ios::sync_with_stdio(false);
  std::cin.tie(nullptr);

  chrom initial = read_file();
  population initial_population =
      generate_init_population(POPULATION_SIZE, initial);

  for (int i = 0; i != MAX_POPULATIONS; ++i) {
  }

  return 0;
}

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

gen generate_gen(gen& initial) {
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

chrom generate_chromosome(chrom& initial) {
  chrom new_chrom;
  for (int row = 0; row != GRID_SIZE; ++row) {
    new_chrom.push_back(generate_gen(initial[row]));
  }
  return new_chrom;
}

population generate_init_population(int size, chrom& initial) {
  population population;
  for (int i = 0; i != size; ++i) {
    population.push_back(generate_chromosome(initial));
  }
  return population;
}

int fitness(chrom& current) {
  int total_fit = 0;

  // Validate rows and columns
  for (int idx = 0; idx != GRID_SIZE; ++idx) {
    int row_sum = 0, col_sum = 0;
    std::vector<bool> undef_row(GRID_SIZE, true);
    std::vector<bool> undef_col(GRID_SIZE, true);

    for (int col = 0; col != GRID_SIZE; ++col) {
      row_sum += current[idx][col].first;

      if (undef_col[current[idx][col].first - 1]) {
        undef_col[current[idx][col].first - 1] = false;
      }
    }
    for (int row = 0; row != GRID_SIZE; ++row) {
      col_sum += current[row][idx].first;

      if (undef_col[current[row][idx].first - 1]) {
        undef_col[current[row][idx].first - 1] = false;
      }
    }

    int row_total_undef =
        std::accumulate(undef_row.begin(), undef_row.end(), 0);
    int col_total_undef =
        std::accumulate(undef_row.begin(), undef_row.end(), 0);

    row_sum -= 45;
    col_sum -= 45;

    total_fit += row_total_undef + (row_total_undef + 1) / 2 + col_total_undef +
                 (col_total_undef + 1) / 2;
    total_fit += row_sum + (row_sum + 1) / 2 + col_sum + (col_sum + 1) / 2;
  }

  // Validate blocks
  for (int block_row = 0; block_row != 3; ++block_row) {
    for (int block_col = 0; block_col != 3; ++block_col) {
      int block_sum = 0;
      for (int row = block_row * 3; row != (block_row + 1) * 3; ++row) {
        for (int col = block_col * 3; col != (block_col + 1) * 3; ++col) {
          block_sum += current[row][col].first;
        }
      }
      block_sum -= 45;
      total_fit += block_sum + (block_sum + 1) / 2;
    }
  }

  return total_fit;
}
