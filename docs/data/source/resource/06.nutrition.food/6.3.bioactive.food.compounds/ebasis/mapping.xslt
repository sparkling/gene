<?xml version="1.0" encoding="UTF-8"?>
<!--
  ============================================================
  eBASIS -> Unified Bioactive Compounds Schema Mapping
  ============================================================
  Source: ./schema.json
  Target: ../schema.json (Unified Bioactive Compounds Schema)
  XSLT Version: 3.0

  TRANSFORMATION SUMMARY:
  ┌─────────────────────────┬─────────────────┬──────────────┐
  │ Source Field            │ Target Field    │ Transform    │
  ├─────────────────────────┼─────────────────┼──────────────┤
  │ food.food_name          │ food_name       │ Direct       │
  │ component.component_name│ compound_name   │ Direct       │
  │ component.component_group│compound_class  │ Direct       │
  │ composition_data.value  │ concentration   │ Direct       │
  │ quality_assessment      │ quality_score   │ Flatten      │
  │ (computed)              │ _source         │ Metadata     │
  └─────────────────────────┴─────────────────┴──────────────┘

  NOTES:
  - eBASIS is bioactive compound focused
  - Quality scores follow EuroFIR methodology
-->
<xsl:stylesheet version="3.0"
    xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
    xmlns:xs="http://www.w3.org/2001/XMLSchema"
    xmlns:fn="http://www.w3.org/2005/xpath-functions"
    xmlns:map="http://www.w3.org/2005/xpath-functions/map"
    xmlns:array="http://www.w3.org/2005/xpath-functions/array"
    exclude-result-prefixes="xs fn map array">

  <xsl:output method="json" indent="yes"/>

  <xsl:template name="xsl:initial-template">
    <xsl:param name="input" as="xs:string"/>
    <xsl:variable name="source" select="json-to-xml($input)"/>
    <xsl:apply-templates select="$source/fn:map"/>
  </xsl:template>

  <xsl:template match="fn:map">
    <fn:map>
      <!-- Record ID -->
      <fn:number key="record_id">
        <xsl:value-of select="fn:number[@key='composition_id']"/>
      </fn:number>

      <!-- Food information -->
      <fn:map key="food">
        <fn:string key="name">
          <xsl:value-of select="fn:map[@key='food']/fn:string[@key='food_name']"/>
        </fn:string>
        <xsl:if test="fn:map[@key='food']/fn:string[@key='scientific_name']">
          <fn:string key="scientific_name">
            <xsl:value-of select="fn:map[@key='food']/fn:string[@key='scientific_name']"/>
          </fn:string>
        </xsl:if>
        <xsl:if test="fn:map[@key='food']/fn:string[@key='food_group']">
          <fn:string key="category">
            <xsl:value-of select="fn:map[@key='food']/fn:string[@key='food_group']"/>
          </fn:string>
        </xsl:if>
      </fn:map>

      <!-- Compound information -->
      <fn:map key="compound">
        <fn:string key="name">
          <xsl:value-of select="fn:map[@key='component']/fn:string[@key='component_name']"/>
        </fn:string>
        <fn:string key="class">
          <xsl:value-of select="fn:map[@key='component']/fn:string[@key='component_group']"/>
        </fn:string>
        <xsl:if test="fn:map[@key='component']/fn:string[@key='component_subgroup']">
          <fn:string key="subclass">
            <xsl:value-of select="fn:map[@key='component']/fn:string[@key='component_subgroup']"/>
          </fn:string>
        </xsl:if>
        <xsl:if test="fn:map[@key='component']/fn:string[@key='cas_number']">
          <fn:string key="cas_number">
            <xsl:value-of select="fn:map[@key='component']/fn:string[@key='cas_number']"/>
          </fn:string>
        </xsl:if>
        <xsl:if test="fn:map[@key='component']/fn:number[@key='pubchem_cid']">
          <fn:number key="pubchem_cid">
            <xsl:value-of select="fn:map[@key='component']/fn:number[@key='pubchem_cid']"/>
          </fn:number>
        </xsl:if>
      </fn:map>

      <!-- Concentration data -->
      <fn:map key="concentration">
        <fn:number key="value">
          <xsl:value-of select="fn:map[@key='composition_data']/fn:number[@key='value']"/>
        </fn:number>
        <fn:string key="unit">
          <xsl:value-of select="fn:map[@key='composition_data']/fn:string[@key='unit']"/>
        </fn:string>
        <xsl:if test="fn:map[@key='composition_data']/fn:number[@key='n_samples']">
          <fn:number key="sample_count">
            <xsl:value-of select="fn:map[@key='composition_data']/fn:number[@key='n_samples']"/>
          </fn:number>
        </xsl:if>
      </fn:map>

      <!-- Quality score -->
      <xsl:if test="fn:map[@key='quality_assessment']/fn:number[@key='overall_score']">
        <fn:number key="quality_score">
          <xsl:value-of select="fn:map[@key='quality_assessment']/fn:number[@key='overall_score']"/>
        </fn:number>
      </xsl:if>

      <!-- Source metadata -->
      <fn:map key="_source">
        <fn:string key="database">eBASIS</fn:string>
        <fn:number key="original_id">
          <xsl:value-of select="fn:number[@key='composition_id']"/>
        </fn:number>
        <xsl:if test="fn:map[@key='reference']/fn:number[@key='pubmed_id']">
          <fn:number key="pubmed_id">
            <xsl:value-of select="fn:map[@key='reference']/fn:number[@key='pubmed_id']"/>
          </fn:number>
        </xsl:if>
      </fn:map>

    </fn:map>
  </xsl:template>

</xsl:stylesheet>
