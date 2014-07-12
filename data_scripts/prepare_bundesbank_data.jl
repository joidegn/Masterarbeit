# This file can be used to pre-clean and extract bundesbank data series from the data/folder given a json file resulting from scraping the bundesbank website


using JSON

load_set(id::String) = readcsv(string("../data/bundesbank/", id, ".csv"))[6:end, :]

function get_set_by_id(id, data_sets)
    filter(set->set["id"]==id, data_sets)
end
get_set_by_id(id) = get_set_by_id(id, json)

function periodicity(data_set)
    data = load_set(data_set["id"])
    return length(split(string(data[1,1]), "-"))  # returns 1 for yearly, 2 for monthly, 3 for daily
end

function longest(data_sets, n=1)
    datas = [load_set(data_set["id"]) for data_set in data_sets]
    data_sets[sortperm(datas, by=set->size(set, 1))[end-(n-1):end]]
end
longest(data_sets) = longest(data_sets, 1)[1]

# only return data sets which have a selected flag in the json file (have to be put in manually)
function filter_data_sets(data_sets::Array{Any, 1}, fixed_periodicity::Int64=2, minimum_length=100)  # periodicity can be 1: yearly, 2:monthly, 3: daily
    selected_ids = Array(String, (0))
    for data_sets in json  # different sets with the same name
        selected = false
        for data_set in data_sets
            if "selected" in keys(data_set)
                selected = true
            end
        end
        if selected
            filtered_set = filter(data_set->periodicity(data_set)==fixed_periodicity && size(load_set(data_set["id"]), 1)>=minimum_length, data_sets)
            if length(filtered_set) > 1
                append!(selected_ids, [longest(filtered_set)["id"]]) # append only the id of the first data set with monthly data
            end
        end
    end
    return(selected_ids)
end



file = open("../data/items.grouped.json")
json = JSON.parse(file)
selected_ids = filter_data_sets(json)

# get the selected series in a seperate folder for easier manual inspection and for the R script used for demeaning: seperate folder would have to be created first
#for id in selected_ids
#    cp(string("data/bundesbank/", id, ".csv"), string("data/bundesbank_selected/", id, ".csv"))
#end

data_sets = [load_set(id) for id in selected_ids]
lengths = [size(data_set, 1) for data_set in data_sets]


length_order = sortperm(data_sets, by=length)
# lengths[length_order]
selected_ids[length_order[1:5]]  # gives the first 5 shortes series
