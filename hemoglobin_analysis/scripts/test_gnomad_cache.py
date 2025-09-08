#!/usr/bin/env python3
"""
Test script for gnomAD caching and backoff functionality
"""

import json
import os
import time
import requests

def test_gnomad_caching():
    """Test gnomAD caching functionality"""
    print("🧪 Testing gnomAD caching and backoff...")

    # Simulate a simple cache
    gnomad_cache = {}

    def query_gnomad_with_caching(variant_id: str) -> dict:
        """Simulate gnomAD query with caching and backoff"""
        cache_key = f"gnomad_{variant_id}"

        # Check cache first
        if cache_key in gnomad_cache:
            cached_result = gnomad_cache[cache_key]
            print(f"📋 Using cached gnomAD result for {variant_id}")
            return cached_result

        # Simulate API call (using a known variant ID for testing)
        url = "https://gnomad.broadinstitute.org/api"
        query = """
        query VariantInfo($variantId: String!) {
            variant(variantId: $variantId) {
                variantId
                exome { af }
                genome { af }
            }
        }
        """

        max_retries = 3
        backoff_delay = 1

        for attempt in range(max_retries):
            try:
                print(f"🌐 Querying gnomAD for {variant_id} (attempt {attempt + 1})")
                response = requests.post(
                    url,
                    json={'query': query, 'variables': {'variantId': variant_id}},
                    timeout=15
                )

                if response.status_code == 200:
                    data = response.json()
                    if data.get('data', {}).get('variant'):
                        variant = data['data']['variant']
                        result = {
                            'found': True,
                            'variant_id': variant_id,
                            'af_exome': variant.get('exome', {}).get('af', 0),
                            'af_genome': variant.get('genome', {}).get('af', 0)
                        }
                        gnomad_cache[cache_key] = result
                        return result
                    else:
                        result = {
                            'found': False,
                            'variant_id': variant_id,
                            'af_exome': 0,
                            'af_genome': 0,
                            'error': 'Variant not found'
                        }
                        gnomad_cache[cache_key] = result
                        return result

                elif response.status_code == 429:
                    delay = backoff_delay * (2 ** attempt)
                    print(f"⏱️ Rate limited, waiting {delay}s before retry")
                    time.sleep(delay)
                    continue

                else:
                    print(f"❌ API error {response.status_code}")
                    if attempt < max_retries - 1:
                        time.sleep(backoff_delay)
                        continue

            except Exception as e:
                print(f"❌ Error: {e}")
                if attempt < max_retries - 1:
                    time.sleep(backoff_delay)
                    continue

        # All retries failed
        result = {
            'found': False,
            'variant_id': variant_id,
            'af_exome': 0,
            'af_genome': 0,
            'error': 'Failed after retries'
        }
        gnomad_cache[cache_key] = result
        return result

    # Test with known hemoglobin variants
    test_variant_ids = [
        "11-5248232-G-A",  # HBB Glu6Val (sickle cell)
        "11-5248232-G-T",  # HBB Glu6Lys (HbC)
    ]

    print(f"\n📊 Initial cache size: {len(gnomad_cache)} entries")

    for variant_id in test_variant_ids:
        print(f"\n🔍 Testing {variant_id}")

        # First query
        result1 = query_gnomad_with_caching(variant_id)
        print(f"First query: found={result1['found']}, exome_af={result1.get('af_exome', 'N/A')}")

        # Second query (should use cache)
        result2 = query_gnomad_with_caching(variant_id)
        print(f"Second query: found={result2['found']}, exome_af={result2.get('af_exome', 'N/A')}")

        if result1 == result2:
            print("✅ Cache working - results identical")
        else:
            print("❌ Cache issue - results differ")

    print(f"\n📊 Final cache size: {len(gnomad_cache)} entries")

    # Save cache to file
    cache_file = "../data/cache/gnomad_cache.json"
    os.makedirs(os.path.dirname(cache_file), exist_ok=True)

    with open(cache_file, 'w') as f:
        json.dump(gnomad_cache, f, indent=2)

    print(f"💾 Saved cache to {cache_file}")
    print("✅ gnomAD caching test completed")

if __name__ == "__main__":
    test_gnomad_caching()
