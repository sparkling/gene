<?xml version="1.0" encoding="UTF-8"?>
<!--
  ============================================================
  Wikidata -> Unified Knowledge Base Schema Mapping
  ============================================================
  Source: ./schema.json
  Target: ../schema.json (Unified Knowledge Base Schema)
  XSLT Version: 3.0

  TRANSFORMATION SUMMARY:
  +-------------------------+-----------------+--------------+
  | Source Field            | Target Field    | Transform    |
  +-------------------------+-----------------+--------------+
  | id                      | id              | Direct       |
  | labels.en.value         | title           | Extract      |
  | descriptions.en.value   | description     | Extract      |
  | claims                  | properties      | Flatten      |
  | sitelinks               | wikipedia_links | Extract      |
  +-------------------------+-----------------+--------------+

  NOTES:
  - English labels used as primary
  - Biomedical properties extracted specifically
-->
<xsl:stylesheet version="3.0"
    xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
    xmlns:xs="http://www.w3.org/2001/XMLSchema"
    xmlns:fn="http://www.w3.org/2005/xpath-functions"
    xmlns:map="http://www.w3.org/2005/xpath-functions/map"
    xmlns:array="http://www.w3.org/2005/xpath-functions/array"
    xmlns:local="http://local.functions"
    exclude-result-prefixes="xs fn map array local">

  <xsl:output method="json" indent="yes"/>

  <xsl:template name="xsl:initial-template">
    <xsl:param name="input" as="xs:string"/>
    <xsl:variable name="source" select="json-to-xml($input)"/>
    <xsl:apply-templates select="$source/fn:map"/>
  </xsl:template>

  <xsl:template match="fn:map">
    <fn:map>
      <!-- Primary identifier -->
      <fn:string key="id">
        <xsl:value-of select="fn:string[@key='id']"/>
      </fn:string>

      <!-- Title from English label -->
      <xsl:if test="fn:map[@key='labels']/fn:map[@key='en']/fn:string[@key='value']">
        <fn:string key="title">
          <xsl:value-of select="fn:map[@key='labels']/fn:map[@key='en']/fn:string[@key='value']"/>
        </fn:string>
      </xsl:if>

      <!-- Description -->
      <xsl:if test="fn:map[@key='descriptions']/fn:map[@key='en']/fn:string[@key='value']">
        <fn:string key="description">
          <xsl:value-of select="fn:map[@key='descriptions']/fn:map[@key='en']/fn:string[@key='value']"/>
        </fn:string>
      </xsl:if>

      <!-- Entity type -->
      <fn:string key="type">
        <xsl:value-of select="fn:string[@key='type']"/>
      </fn:string>

      <!-- Cross-references from claims -->
      <fn:map key="cross_references">
        <!-- NCBI Gene ID (P351) -->
        <xsl:if test="fn:map[@key='claims']/fn:array[@key='P351']">
          <fn:string key="ncbi_gene">
            <xsl:value-of select="fn:map[@key='claims']/fn:array[@key='P351']/fn:map[1]/fn:map[@key='mainsnak']/fn:map[@key='datavalue']/fn:string[@key='value']"/>
          </fn:string>
        </xsl:if>
        <!-- UniProt (P352) -->
        <xsl:if test="fn:map[@key='claims']/fn:array[@key='P352']">
          <fn:string key="uniprot">
            <xsl:value-of select="fn:map[@key='claims']/fn:array[@key='P352']/fn:map[1]/fn:map[@key='mainsnak']/fn:map[@key='datavalue']/fn:string[@key='value']"/>
          </fn:string>
        </xsl:if>
        <!-- MeSH (P486) -->
        <xsl:if test="fn:map[@key='claims']/fn:array[@key='P486']">
          <fn:string key="mesh">
            <xsl:value-of select="fn:map[@key='claims']/fn:array[@key='P486']/fn:map[1]/fn:map[@key='mainsnak']/fn:map[@key='datavalue']/fn:string[@key='value']"/>
          </fn:string>
        </xsl:if>
        <!-- HGNC symbol (P353) -->
        <xsl:if test="fn:map[@key='claims']/fn:array[@key='P353']">
          <fn:string key="hgnc_symbol">
            <xsl:value-of select="fn:map[@key='claims']/fn:array[@key='P353']/fn:map[1]/fn:map[@key='mainsnak']/fn:map[@key='datavalue']/fn:string[@key='value']"/>
          </fn:string>
        </xsl:if>
      </fn:map>

      <!-- Wikipedia links -->
      <xsl:if test="fn:map[@key='sitelinks']/fn:map[@key='enwiki']">
        <fn:string key="wikipedia_url">
          <xsl:value-of select="concat('https://en.wikipedia.org/wiki/', encode-for-uri(fn:map[@key='sitelinks']/fn:map[@key='enwiki']/fn:string[@key='title']))"/>
        </fn:string>
      </xsl:if>

      <!-- Source metadata -->
      <fn:map key="_source">
        <fn:string key="database">Wikidata</fn:string>
        <fn:string key="original_id">
          <xsl:value-of select="fn:string[@key='id']"/>
        </fn:string>
        <fn:string key="url">
          <xsl:value-of select="concat('https://www.wikidata.org/wiki/', fn:string[@key='id'])"/>
        </fn:string>
      </fn:map>

    </fn:map>
  </xsl:template>

</xsl:stylesheet>
