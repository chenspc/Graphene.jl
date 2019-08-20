export nn2dfb

# function nn2dfb(nn_path::String, graphene_path::String, dfb_path::String; dataset=[], range=[])
function nn2dfb(nn_path::String, dfb_path::String; dataset=[], range=[])
    if dataset == []
        dataset = nn_path
    end
    xy_stack = import_stack(nn_path; range=range)
    graphene_stack = map(x -> make_graphene(xy_stack[x]; frame=x, dataset=dataset), 1:length(xy_stack))
    # FileIO.save(graphene_path, merge_stack(graphene_stack))
    find_dfb_save(graphene_stack, dfb_path)
    # return 1
end
