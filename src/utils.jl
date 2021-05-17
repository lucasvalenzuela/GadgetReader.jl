skipnan(x) = Iterators.filter(!isnan, x)

function get_dims(perspective::Symbol)
    dims = if perspective === :edgeon
        [1, 3]
    elseif perspective === :sideon
        [2, 3]
    elseif perspective === :faceon
        [1, 2]
    else
        error("`perspective` is $perspective, needs to be one of `:edgeon`, `:sideon`, `:faceon`")
    end
end
