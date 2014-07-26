using Fred
using JSON

const API_KEY = "e3a3139530d8a6ef68e36b8b12bbc2f7"
const BASE_URL = "http://api.stlouisfed.org/fred/"

ids = readcsv("../data/series_FRED")
ids = convert(Array{ASCIIString, 1}, reshape(ids, size(ids, 1)))

get_request(uri, query) = get("$(BASE_URL)$(uri)"; query=Base.merge({"api_key"=>API_KEY, "file_type"=>"json"}, query))  # request with API_KEY and file_type
requests = [get_request("series", {"series_id"=>series_id}) for series_id in ids]
jsons = [JSON.parse(req.data) for req in requests]
jsons_success = jsons[Bool["seriess" in keys(json) for json in jsons]]
jsons_failure = jsons[Bool[!in("seriess", keys(json)) for json in jsons]]  # should be 0
titles = [json["seriess"][1]["title"] for json in jsons]

#file = open("../data/series_FRED.tex", "w+")
#write(file, matrix_to_table(hcat(ids, titles)))
#close(file)

#writecsv("series_FRED_descriptions", hcat(ids, titles))
