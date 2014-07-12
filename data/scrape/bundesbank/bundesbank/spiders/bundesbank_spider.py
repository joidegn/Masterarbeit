import scrapy

from bundesbank.items import DataSet

class BundesbankSpider(scrapy.Spider):
    name = "bundesbank"
    allowed_domains = ["bundesbank.de"]
    base_domain = "http://www.bundesbank.de"
    start_urls = [
        "http://www.bundesbank.de/Navigation/EN/Statistics/Time_series_databases/Macro_economic_time_series/macro_economic_time_series_node.html?openAll=true"
    ]

    def parse(self, response):
        for url in response.selector.xpath("//a[@class='list']/@href").extract():
            yield scrapy.Request(url, callback=self.parse_csv_list)

    def parse_csv_list(self, response):
        for tr in response.selector.xpath("//tr"):
            urls = tr.xpath("td/a[contains(@href, 'its_csvFormat=en&its_fileFormat=csv')]/@href").extract()
            if len(urls) > 0:
                url = self.base_domain + "/" + urls[0]
                data_set = DataSet()
                data_set["id"] = tr.xpath("td[1]/a/text()").extract()[0].strip()
                data_set["title"] = tr.xpath("//td[2]/text()").extract()[0].strip()
                yield scrapy.Request(url, meta={"data_set": data_set}, callback=self.parse_csv)


    def parse_csv(self, response):
        data_set = response.request.meta['data_set']
        data_set["csv_url"] = response.url
        with open("data/%s.csv" % data_set["id"], "wb") as f:
            f.write(response.body)
        yield data_set
