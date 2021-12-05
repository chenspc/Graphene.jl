using Graphene
using Test

# make_gsa
# link_past_future!
# chain_time_series!
# plot_graphene_gsa
# save2gif

test_xy1 = [0     0     8.66    8.66    17.32    17.32;
            5     15    0       20      5        15]
test_graphene1 = make_graphene(test_xy1; max_bondlength=15, frame=1)


test_shifts = [0.012     0.013    0.014     0.021     0.031     0.041;
               0.021     0.031    0.041     0.012     0.013     0.014]
test_xy2 = test_xy1 + test_shifts
test_graphene2 = make_graphene(test_xy2; max_bondlength=15, frame=2)

test_gsa1 = make_gsa(test_graphene1)
test_gsa2 = make_gsa(test_graphene2)

link_past_future!(test_gsa1, test_gsa2)
# find_past_future!(test_graphene1, test_graphene2)
test_gsa1
test_gsa2

t

size(test_gsa1)
size(t)
