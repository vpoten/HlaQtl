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
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import groovy.json.JsonOutput;
import org.json.JSONArray;
import org.json.JSONObject;

import org.ngsutils.variation.SNPData;

/**
 *
 * @author victor
 */
public class EnsemblRestClient {
    private final String server;
    private final int requestPerSecond;
    private final int maxPost;
    private int requestCount = 0;
    private long limitStartTime = System.currentTimeMillis();

    public EnsemblRestClient(String server, int requestPerSecond, int maxPost) {
        this.server = server;
        this.requestPerSecond = requestPerSecond;
        this.maxPost = maxPost;
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

        EnsemblRestClient client = new EnsemblRestClient("http://rest.ensembl.org", 15, 200);
        client.printVariants(species, symbol);
    }

    public void printVariants(String species, String symbol) throws UnirestException, InterruptedException {
        String geneId = getGeneId(species, symbol);

        String url = String.format("%s/overlap/id/%s?feature=variation", server, geneId);
        JSONArray variants = fetchJson(url, "get", null).getArray();
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
        JSONArray genes = fetchJson(url, "get", null).getArray();

        if (genes.length() == 0) {
            throw new IllegalArgumentException(String.format("No gene id for %s symbol %s", species, symbol));
        } else {
            return genes.getJSONObject(0).getString("id");
        }
    }
    
    private SNPData parseJSONObjectSNP(JSONObject snpObj) {
        SNPData snp = new SNPData();
        snp.setId(snpObj.getString("name"));
        snp.setMaf(snpObj.getDouble("MAF"));
        snp.setMinor(snpObj.getString("minor_allele"));
        
        JSONObject mapping = snpObj.getJSONArray("mappings").getJSONObject(0);
        
        snp.setChr(mapping.getString("seq_region_name"));
        snp.setPosition(mapping.getInt("start"));
        snp.setAlleles(mapping.getString("allele_string"));
        
        return snp;
    }
    
    public List<SNPData> getSnps(List<String> ids, String species) throws UnirestException, InterruptedException {
        String url = String.format("%s/variation/%s", server, species);
        
        List result = new ArrayList();
        int start = 0;

        while(start < ids.size()) {
            HashMap<String, List<String>> requestBody = new HashMap<String, List<String>>();
            int end = (start + this.maxPost > ids.size()) ? ids.size() : start + this.maxPost;
            requestBody.put("ids", ids.subList(start, end));
            JSONObject snps = fetchJson(url, "post", JsonOutput.toJson(requestBody)).getObject();
            
            for(Object key : snps.keySet()) {
                result.add(parseJSONObjectSNP(snps.getJSONObject((String)key)));
            }
            
            start +=  snps.length();
        }

        return result;
    }

    private JsonNode fetchJson(String url, String method, String body) throws UnirestException, InterruptedException {
        rateLimit();
        HttpResponse<JsonNode> response;
        
        if (method.compareTo("get") == 0) {
            response = Unirest.get(url)
                    .header("Content-Type", "application/json")
                    .asJson();
        } else {
            response = Unirest.post(url)
                    .header("Accept", "application/json")
                    .header("Content-Type", "application/json")
                    .body(body)
                    .asJson();
        }
        
        String retryHeader = response.getHeaders().getFirst("Retry-After");

        if (response.getStatus() == 200) {
            return response.getBody();
        } else if (response.getStatus() == 429 && retryHeader != null) {
            Long waitSeconds = Long.valueOf(retryHeader);
            Thread.sleep(waitSeconds * 1000);
            return fetchJson(url, method, body);
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
