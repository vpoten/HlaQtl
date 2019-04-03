/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.msgenetics.hlaqtl;

import com.mashape.unirest.http.HttpResponse;
import com.mashape.unirest.http.JsonNode;
import com.mashape.unirest.http.Unirest;
import com.mashape.unirest.http.exceptions.UnirestException;
import org.json.JSONArray;
import org.json.JSONObject;


/**
 *
 * @author victor
 */
public class EnsemblRestClient {
    private final String server;
    private final int requestPerSecond;
    private int requestCount = 0;
    private long limitStartTime = System.currentTimeMillis();

    public EnsemblRestClient(String server, int requestPerSecond) {
        this.server = server;
        this.requestPerSecond = requestPerSecond;
    }

    public static void main(String[] args) throws InterruptedException, UnirestException {
        String species = "human";
        String symbol = "BRAF";
        if (args.length == 1) {
            species = args[0];
        } else if (args.length == 2) {
            species = args[0];
            symbol = args[1];
        }

        EnsemblRestClient client = new EnsemblRestClient("http://rest.ensembl.org", 15);
        client.printVariants(species, symbol);
    }

    public void printVariants(String species, String symbol) throws UnirestException, InterruptedException {
        String geneId = getGeneId(species, symbol);

        String url = String.format("%s/overlap/id/%s?feature=variation", server, geneId);
        JSONArray variants = fetchJson(url).getArray();
        for (int i = 0; i < variants.length(); i++) {
            JSONObject variant = variants.getJSONObject(i);
            String srName = variant.getString("seq_region_name");
            int start = variant.getInt("start");
            int end = variant.getInt("end");
            int strand = variant.getInt("strand");
            String id = variant.getString("id");
            String consequence = variant.getString("consequence_type");
            String output = String.format("%s:%d-%d:%d ==> %s (%s)", srName, start, end, strand, id, consequence);
            System.out.println(output);
        }
    }

    private String getGeneId(String species, String symbol) throws UnirestException, InterruptedException {
        String url = String.format("%s/xrefs/symbol/%s/%s?object_type=gene", server, species, symbol);
        JSONArray genes = fetchJson(url).getArray();

        if (genes.length() == 0) {
            throw new IllegalArgumentException(String.format("No gene id for %s symbol %s", species, symbol));
        } else {
            return genes.getJSONObject(0).getString("id");
        }
    }

    private JsonNode fetchJson(String url) throws UnirestException, InterruptedException {
        rateLimit();
        HttpResponse<JsonNode> response = Unirest.get(url)
                        .header("Content-Type", "application/json")
                        .asJson();
        String retryHeader = response.getHeaders().getFirst("Retry-After");

        if (response.getStatus() == 200) {
            return response.getBody();
        } else if (response.getStatus() == 429 && retryHeader != null) {
            Long waitSeconds = Long.valueOf(retryHeader);
            Thread.sleep(waitSeconds * 1000);
            return fetchJson(url);
        } else {
            throw new IllegalArgumentException("No data at " + url);
        }
    }

    private void rateLimit() throws InterruptedException {
        requestCount++;
        if (requestCount == requestPerSecond) {
            long currentTime = System.currentTimeMillis();
            long diff = currentTime - limitStartTime;
            //if less than a second has passed then sleep for the remainder of the second
            if (diff < 1000) {
                    Thread.sleep(1000 - diff);
            }
            //reset
            limitStartTime = System.currentTimeMillis();
            requestCount = 0;
        }
    }
}
