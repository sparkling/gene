<?xml version="1.0" encoding="UTF-8"?>
<!--
  ============================================================
  Open Food Facts -> Unified Food Composition Schema Mapping
  ============================================================
  Source: ./schema.json
  Target: ../schema.json (Unified Food Composition Schema)
  XSLT Version: 3.0

  TRANSFORMATION SUMMARY:
  ┌─────────────────────────┬─────────────────┬──────────────┐
  │ Source Field            │ Target Field    │ Transform    │
  ├─────────────────────────┼─────────────────┼──────────────┤
  │ code                    │ food_id         │ Direct       │
  │ product_name            │ food_name       │ Direct       │
  │ brands                  │ brand           │ Split array  │
  │ categories_tags         │ food_category   │ Direct       │
  │ nutriments              │ nutrient_values │ Restructure  │
  │ nutriscore_grade        │ nutriscore_grade│ Direct       │
  │ nova_group              │ nova_group      │ Direct       │
  │ ecoscore_grade          │ ecoscore_grade  │ Direct       │
  │ quantity                │ quantity        │ Direct       │
  │ (computed)              │ _source         │ Metadata     │
  └─────────────────────────┴─────────────────┴──────────────┘

  NULL HANDLING:
  - product_name: Use generic_name if empty
  - brands: Split on comma if present

  NOTES:
  - Open Food Facts uses barcode as primary identifier
  - Nutrient values are per 100g/100ml
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

  <!-- ============================================================
       Entry Point: Parse JSON input
       ============================================================ -->
  <xsl:template name="xsl:initial-template">
    <xsl:param name="input" as="xs:string"/>
    <xsl:variable name="source" select="json-to-xml($input)"/>
    <xsl:apply-templates select="$source/fn:map"/>
  </xsl:template>

  <!-- ============================================================
       Main Record Transformation
       ============================================================ -->
  <xsl:template match="fn:map">
    <fn:map>
      <!-- ========== DIRECT FIELD MAPPINGS ========== -->

      <!-- Primary identifier: barcode -->
      <fn:string key="food_id">
        <xsl:value-of select="fn:string[@key='code']"/>
      </fn:string>

      <!-- Product name -->
      <fn:string key="food_name">
        <xsl:choose>
          <xsl:when test="fn:string[@key='product_name'] != ''">
            <xsl:value-of select="fn:string[@key='product_name']"/>
          </xsl:when>
          <xsl:when test="fn:string[@key='generic_name'] != ''">
            <xsl:value-of select="fn:string[@key='generic_name']"/>
          </xsl:when>
          <xsl:otherwise>Unknown Product</xsl:otherwise>
        </xsl:choose>
      </fn:string>

      <!-- Scientific name: not applicable for packaged products -->
      <fn:null key="scientific_name"/>

      <!-- Food category: pass through categories_tags -->
      <xsl:choose>
        <xsl:when test="fn:array[@key='categories_tags']">
          <fn:array key="food_category">
            <xsl:for-each select="fn:array[@key='categories_tags']/fn:string">
              <fn:string><xsl:value-of select="."/></fn:string>
            </xsl:for-each>
          </fn:array>
        </xsl:when>
        <xsl:otherwise>
          <fn:array key="food_category"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- Brand: split comma-separated into array -->
      <xsl:choose>
        <xsl:when test="fn:string[@key='brands'] and fn:string[@key='brands'] != ''">
          <fn:array key="brand">
            <xsl:for-each select="tokenize(fn:string[@key='brands'], ',')">
              <fn:string><xsl:value-of select="normalize-space(.)"/></fn:string>
            </xsl:for-each>
          </fn:array>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="brand"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== NUTRIENT VALUES ========== -->

      <xsl:if test="fn:map[@key='nutriments']">
        <fn:map key="nutrient_values">
          <xsl:variable name="nutr" select="fn:map[@key='nutriments']"/>

          <xsl:if test="$nutr/fn:number[@key='energy-kcal_100g']">
            <fn:number key="energy_kcal">
              <xsl:value-of select="$nutr/fn:number[@key='energy-kcal_100g']"/>
            </fn:number>
          </xsl:if>

          <xsl:if test="$nutr/fn:number[@key='proteins_100g']">
            <fn:number key="protein_g">
              <xsl:value-of select="$nutr/fn:number[@key='proteins_100g']"/>
            </fn:number>
          </xsl:if>

          <xsl:if test="$nutr/fn:number[@key='fat_100g']">
            <fn:number key="fat_g">
              <xsl:value-of select="$nutr/fn:number[@key='fat_100g']"/>
            </fn:number>
          </xsl:if>

          <xsl:if test="$nutr/fn:number[@key='carbohydrates_100g']">
            <fn:number key="carbs_g">
              <xsl:value-of select="$nutr/fn:number[@key='carbohydrates_100g']"/>
            </fn:number>
          </xsl:if>

          <xsl:if test="$nutr/fn:number[@key='fiber_100g']">
            <fn:number key="fiber_g">
              <xsl:value-of select="$nutr/fn:number[@key='fiber_100g']"/>
            </fn:number>
          </xsl:if>

          <xsl:if test="$nutr/fn:number[@key='sugars_100g']">
            <fn:number key="sugar_g">
              <xsl:value-of select="$nutr/fn:number[@key='sugars_100g']"/>
            </fn:number>
          </xsl:if>

          <xsl:if test="$nutr/fn:number[@key='salt_100g']">
            <fn:number key="sodium_mg">
              <xsl:value-of select="$nutr/fn:number[@key='salt_100g'] * 400"/>
            </fn:number>
          </xsl:if>
        </fn:map>
      </xsl:if>

      <!-- ========== SCORES ========== -->

      <!-- Nutri-Score -->
      <xsl:choose>
        <xsl:when test="fn:string[@key='nutriscore_grade']">
          <fn:string key="nutriscore_grade">
            <xsl:value-of select="fn:string[@key='nutriscore_grade']"/>
          </fn:string>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="nutriscore_grade"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- NOVA group -->
      <xsl:choose>
        <xsl:when test="fn:number[@key='nova_group']">
          <fn:number key="nova_group">
            <xsl:value-of select="fn:number[@key='nova_group']"/>
          </fn:number>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="nova_group"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- Eco-Score -->
      <xsl:choose>
        <xsl:when test="fn:string[@key='ecoscore_grade']">
          <fn:string key="ecoscore_grade">
            <xsl:value-of select="fn:string[@key='ecoscore_grade']"/>
          </fn:string>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="ecoscore_grade"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- Quantity -->
      <xsl:choose>
        <xsl:when test="fn:string[@key='quantity']">
          <fn:string key="quantity">
            <xsl:value-of select="fn:string[@key='quantity']"/>
          </fn:string>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="quantity"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- Not applicable for OFF -->
      <fn:null key="compound_content"/>
      <fn:null key="moldb_smiles"/>
      <fn:null key="moldb_inchikey"/>
      <fn:null key="cas_number"/>
      <fn:null key="itis_id"/>

      <!-- ========== METADATA ========== -->

      <fn:map key="_source">
        <fn:string key="database">Open Food Facts</fn:string>
        <fn:string key="version">2026-01</fn:string>
        <fn:string key="original_id">
          <xsl:value-of select="fn:string[@key='code']"/>
        </fn:string>
      </fn:map>

    </fn:map>
  </xsl:template>

</xsl:stylesheet>
