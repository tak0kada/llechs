all: calc_coef recv_shape post_acvd

calc_coef:
	g++ -std=c++17 -Wall -Wextra -Werror -Wno-unused-value -Wno-unused-local-typedefs -Wno-long-long -Wno-variadic-macros -fexceptions -O2 -I. -I./include main_calc_coef.cpp -lboost_program_options -o calc_coef

recv_shape:
	g++ -std=c++17 -Wall -Wextra -Werror -Wno-unused-value -Wno-unused-local-typedefs -Wno-long-long -Wno-variadic-macros -fexceptions -O2 -I. -I./include main_recv_shape.cpp -lboost_program_options -o recv_shape

post_acvd:
	g++ -std=c++17 -Wall -Wextra -Werror -Wno-unused-value -Wno-unused-local-typedefs -Wno-long-long -Wno-variadic-macros -fexceptions -O2 -I. -I./include main_post_acvd.cpp -lboost_program_options -o post_acvd
