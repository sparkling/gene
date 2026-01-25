<?xml version="1.0" encoding="UTF-8"?>
<!--
  ============================================================
  USDA FoodData Central -> Food Compounds Unified Schema Mapping
  ============================================================
  Source: ./schema.md (USDA FoodData Central)
  Target: ../schema.json (2.4 Food Compounds and Nutrients Unified Schema)
  XSLT Version: 3.0

  TRANSFORMATION SUMMARY:
  ┌─────────────────────────────┬─────────────────────┬──────────────┐
  │ Source Field                │ Target Field        │ Transform    │
  ├─────────────────────────────┼─────────────────────┼──────────────┤
  │ fdcId                       │ fdc_id              │ Direct       │
  │ fdcId                       │ compound_id         │ String       │
  │ description                 │ name                │ Direct       │
  │ dataType                    │ data_type           │ Direct       │
  │ nutrient_id                 │ nutrient_id         │ Direct       │
  │ nutrientNumber              │ nutrient_nbr        │ Direct       │
  │ foodPortions                │ food_portions       │ Object array │
  │ gtinUpc                     │ gtin_upc            │ Direct       │
  │ amount                      │ content_value       │ Direct       │
  │ unitName                    │ content_unit        │ Direct       │
  └─────────────────────────────┴─────────────────────┴──────────────┘

  NULL HANDLING:
  - description: "" → null
  - All numeric fields: null/undefined → null

  NOTES:
  - FDC ID is the primary identifier
  - Data types: Foundation, SR Legacy, Branded, FNDDS
  - Includes serving size conversions
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
      <!-- ========== IDENTIFIERS ========== -->

      <xsl:choose>
        <xsl:when test="fn:number[@key='fdcId']">
          <fn:number key="fdc_id">
            <xsl:value-of select="fn:number[@key='fdcId']"/>
          </fn:number>
          <fn:string key="compound_id">
            <xsl:value-of select="fn:number[@key='fdcId']"/>
          </fn:string>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="fdc_id"/>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="fn:number[@key='nutrient_id']">
          <fn:number key="nutrient_id">
            <xsl:value-of select="fn:number[@key='nutrient_id']"/>
          </fn:number>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="nutrient_id"/>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='nutrientNumber'])">
          <fn:null key="nutrient_nbr"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="nutrient_nbr">
            <xsl:value-of select="fn:string[@key='nutrientNumber']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='gtinUpc'])">
          <fn:null key="gtin_upc"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="gtin_upc">
            <xsl:value-of select="fn:string[@key='gtinUpc']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== BASIC PROPERTIES ========== -->

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='description'])">
          <fn:null key="name"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="name">
            <xsl:value-of select="fn:string[@key='description']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='dataType'])">
          <fn:null key="data_type"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="data_type">
            <xsl:value-of select="fn:string[@key='dataType']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== NUTRIENT CONTENT ========== -->

      <xsl:choose>
        <xsl:when test="fn:number[@key='amount']">
          <fn:number key="content_value">
            <xsl:value-of select="fn:number[@key='amount']"/>
          </fn:number>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="content_value"/>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='unitName'])">
          <fn:null key="content_unit"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="content_unit">
            <xsl:value-of select="fn:string[@key='unitName']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== FOOD PORTIONS ========== -->

      <xsl:choose>
        <xsl:when test="fn:array[@key='foodPortions']/fn:map">
          <fn:array key="food_portions">
            <xsl:for-each select="fn:array[@key='foodPortions']/fn:map">
              <fn:map>
                <xsl:if test="fn:string[@key='measureUnit'] or fn:string[@key='modifier']">
                  <fn:string key="measure_unit">
                    <xsl:value-of select="concat(fn:string[@key='modifier'], ' ', fn:string[@key='measureUnit'])"/>
                  </fn:string>
                </xsl:if>
                <xsl:if test="fn:number[@key='gramWeight']">
                  <fn:number key="gram_weight">
                    <xsl:value-of select="fn:number[@key='gramWeight']"/>
                  </fn:number>
                </xsl:if>
                <xsl:if test="fn:number[@key='amount']">
                  <fn:number key="amount">
                    <xsl:value-of select="fn:number[@key='amount']"/>
                  </fn:number>
                </xsl:if>
              </fn:map>
            </xsl:for-each>
          </fn:array>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="food_portions"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== METADATA ========== -->

      <fn:map key="_source">
        <fn:string key="database">USDA FoodData</fn:string>
        <fn:string key="version">
          <xsl:value-of select="fn:string[@key='source_version']"/>
        </fn:string>
        <fn:string key="original_id">
          <xsl:value-of select="fn:number[@key='fdcId']"/>
        </fn:string>
      </fn:map>

    </fn:map>
  </xsl:template>

  <xsl:function name="local:is-null" as="xs:boolean">
    <xsl:param name="value" as="xs:string?"/>
    <xsl:sequence select="not($value) or $value = '' or $value = '-' or $value = 'NA' or $value = 'N/A' or $value = 'null'"/>
  </xsl:function>

</xsl:stylesheet>
